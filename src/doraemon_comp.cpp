#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

enum h_tag {
  code = 0, lz_ofs = 1, bits = 2
};

struct lha_config {
  size_t header_size;
  size_t command_limit;
  bool use_old_encoding;
  bool allow_empty;
  size_t lz_ofs_bits;
  std::array<size_t, 3> max_bitlen;
  std::function<bool(int, size_t)> allow_zero_bits;
};

std::vector<uint8_t> doraemon_comp_core(
    std::span<const uint8_t> input, const lha_config& config) {

  enum method { uncomp, lz };
  using tag = tag_o<method>;

  struct encodes {
    encode::huffman_t code;
    encode::huffman_t ofs;
    encode::huffman_t bits;
  };

  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = lz_min_len + 0xfd;

  static constexpr size_t huff_offset = 0x0100 - lz_min_len;
  static constexpr size_t lz_ofs_max_bits = 0x0f;
  assert(config.lz_ofs_bits <= lz_ofs_max_bits);

  assert(1 <= config.command_limit && config.command_limit <= 0xffff);
  const size_t chunk_size = config.command_limit;

  const auto encode_codes = [&config](std::span<const size_t> counts, h_tag tag) {
    auto ret = encode::length_limited_huffman(counts, config.max_bitlen[tag], true);
    if (ret.words.size() == 1) {
      if (!(config.use_old_encoding && config.allow_zero_bits(tag, ret.words[0]))) {
        ret.code[ret.words[0]].bitlen += 1;
      }
    }
    return ret;
  };

  const auto max_elem = [](std::span<const size_t> words) -> ptrdiff_t {
    if (words.size() == 0) return -1; // This value must be -1.
    return *std::max_element(words.begin(), words.end());
  };

  const auto bits_table = [&max_elem](const encode::huffman_t& huff) {
    const size_t table_size = max_elem(huff.words) + 1;
    std::vector<uint8_t> bits(table_size, 0);
    for (const auto w : huff.words) bits[w] = huff.code[w].bitlen;
    return bits;
  };

  const auto encode_table = [&bits_table, &encode_codes](const encode::huffman_t& huff) {
    const auto bits = bits_table(huff);
    if (bits.empty()) return encode::huffman_t();

    const size_t max_bit_size = huff.code[huff.words.back()].bitlen;
    std::vector<size_t> counts(3 + max_bit_size, 0);
    for (size_t i = 0; i < bits.size(); ) {
      if (bits[i] != 0) {
        counts[bits[i] + 2] += 1;
        i += 1;
      } else {
        size_t l = encode::run_length(bits, i, 0);
        i += l;
        // greedy (i.e. might not be optimal)
        while (l >= 0x14) counts[2] += 1, l -= std::min<size_t>(l, 0x14 + (1 << 9) - 1);
        while (l >= 0x03) counts[1] += 1, l -= std::min<size_t>(l, 0x03 + (1 << 4) - 1);
        while (l >= 0x01) counts[0] += 1, l -= std::min<size_t>(l, 0x01);
      }
    }
    return encode_codes(counts, h_tag::bits);
  };

  const auto lz_ofs = [&config] {
    size_t k = config.lz_ofs_bits;
    std::vector<vrange> ret(k + 1);
    ret[0] = vrange(1, 1, 0, 0);
    for (size_t k = 1; k < ret.size(); ++k) {
      ret[k] = vrange((1 << (k - 1)) + 1, 1 << k, (k > 0) ? k - 1 : 0, 0);
    }
    return ret;
  }();

  static constexpr auto lz_lens = create_array<size_t, lz_max_len - lz_min_len + 1>([&](size_t i) {
    return i + lz_min_len;
  });

  const auto lz_memo = [&] {
    std::vector<std::array<encode::lz_data, lz_ofs_max_bits + 1>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, lz_ofs, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find(i, max_ofs, lz_min_len);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  using namespace data_type;
  writer_b8_h ret(config.header_size);

  for (size_t begin = 0; begin < input.size(); ) {
    const size_t end = std::min<size_t>(begin + chunk_size, input.size());
    const size_t size = end - begin;

    auto code_cost = [&] {
      // [Todo] Find a reasonable initialization.
      std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
      const auto h = encode::huffman(utility::freq_u8(std::span(input).subspan(begin, size)), true);
      for (const auto w : h.words) ret[w] = h.code[w].bitlen;
      return ret;
    }();
    std::vector<vrange> curr_ofs(lz_ofs.begin(), lz_ofs.end());
    for (auto& v : curr_ofs) v.bitlen += 4;

    using node_type = solver<tag>::node;
    size_t best_cost = std::numeric_limits<size_t>::max();
    std::vector<node_type> best_commands;
    encodes best_enc;

    for (bool updated = true; updated; ) {
      solver<tag> dp(size);

      const auto shift_lz = [&](const encode::lz_data& p) -> encode::lz_data {
        return {p.ofs - begin, p.len}; // Note: can be negative
      };

      for (size_t i = size; i-- > 0; ) {
        dp.update_matrix(i, curr_ofs, lz_lens,
          [&](size_t li) { return code_cost[lz_lens[li] + huff_offset]; },
          [&](size_t oi) { return shift_lz(lz_memo[begin + i][oi]); },
          [&](size_t oi, size_t) -> tag { return {lz, oi}; }
        );
        dp.update(i, 1, code_cost[input[begin + i]], {uncomp, 0});
      }

      struct counter_t {
        std::array<size_t, 0x100 + lz_lens.size()> code = {};
        std::array<size_t, lz_ofs_max_bits + 1> lz_ofs = {};
      } counter;

      const auto commands = [&] {
        std::vector<node_type> ret;
        size_t adr = 0;
        while (adr < size) {
          node_type cmd = dp[adr];
          const auto [tag, oi] = cmd.type;
          if (tag == uncomp) counter.code[input[begin + adr]] += 1;
          else counter.code[cmd.len + huff_offset] += 1, counter.lz_ofs[oi] += 1;
          adr += cmd.len;
          ret.emplace_back(cmd);
        }
        assert(adr == size);
        return ret;
      }();

      const auto calc_cost = [&](const encodes& enc, const counter_t& counter) {
        size_t ret = 0;
        for (const auto& w : enc.code.words) ret += counter.code[w] * enc.code.code[w].bitlen;
        for (const auto& w : enc.ofs.words) ret += counter.lz_ofs[w] * (enc.ofs.code[w].bitlen + lz_ofs[w].bitlen);
        return ret;
      };

      const auto update_costs = [&](const encodes& enc) {
        const auto penalty = [&](const encode::huffman_t& huff, size_t penalty = 0) {
          size_t longest_bitlen = 0;
          for (const auto w : huff.words) longest_bitlen = std::max<size_t>(longest_bitlen, huff.code[w].bitlen);
          return longest_bitlen + penalty;
        };
        code_cost.assign(code_cost.size(), penalty(enc.code, 1));
        for (const auto w : enc.code.words) code_cost[w] = enc.code.code[w].bitlen;

        const size_t c = penalty(enc.ofs);
        for (auto& v : curr_ofs) v.bitlen = c;
        for (const auto w : enc.ofs.words) curr_ofs[w].bitlen = enc.ofs.code[w].bitlen + lz_ofs[w].bitlen;
      };

      const encodes enc {
        .code = encode_codes(counter.code, h_tag::code),
        .ofs = encode_codes(counter.lz_ofs, h_tag::lz_ofs),
        .bits = encode_table(enc.code)
      };

      const size_t estimated_cost = calc_cost(enc, counter);
      if (estimated_cost < best_cost) {
        best_cost = estimated_cost;
        best_commands = std::move(commands);
        best_enc = std::move(enc);
        update_costs(best_enc);
      } else {
        updated = false;
      }
    }

    const auto write_huff_bits = [&](const encode::huffman_t& huff, const size_t s, const size_t t, const size_t lim, const h_tag tag) {
      auto bits = bits_table(huff);
      if (bits.size() > lim) std::runtime_error("This algorithm cannot compress the given data.");
      if (!config.allow_empty) {
        if (huff.words.size() == 0) {
          ret.write<bnh>({s, 1});
          ret.write<bnh>({3, 0});
          return;
        }
      }

      if (config.use_old_encoding) {
        if (huff.words.size() == 0) {
          ret.write<bnh>({s, 1});
          ret.write<bnh>({3, 0});
          return;
        } else if (huff.words.size() == 1) {
          if (config.allow_zero_bits(tag, huff.words[0])) {
            // Caution: These codes might not work on the original LHA algorithm.
            ret.write<bnh>({s, 0});
            ret.write<bnh>({s, huff.words[0]});
            return;
          }
        }
      }

      ret.write<bnh>({s, bits.size()});
      for (size_t i = 0; i < bits.size(); ) {
        if (bits[i] < 7) {
          ret.write<bnh>({3, bits[i]});
        } else {
          ret.write<bnh>({3, 7});
          for (size_t j = bits[i]; j >= 8; --j) ret.write<bnh>({1, 1});
          ret.write<bnh>({1, 0});
        }
        i += 1;
        if (i == t) {
          size_t l = encode::run_length(bits, i, 0);
          l = std::min<size_t>(l, 3);
          if (bits[i] != 0) l = 0;
          ret.write<bnh>({2, l});
          i += l;
        }
      }
    };

    const auto write_table = [&](const encodes& enc) {
      const auto table = bits_table(enc.code);

      const auto& words = enc.code.words;
      if (words.size() == 0) throw std::logic_error("Nothing to compress.");

      if (config.use_old_encoding) {
        if (words.size() == 1) {
          if (config.allow_zero_bits(h_tag::code, words[0])) {
            // Caution: These codes might not work on the original LHA algorithm.
            ret.write<bnh>({9, 0});
            ret.write<bnh>({9, words[0]});
            return;
          }
        }
      }

      const auto write_zeros = [&](size_t l, const size_t z_code, const size_t b, const size_t min_len) -> size_t {
        const size_t max_len = min_len + (1 << b) - 1;
        while (l >= min_len) {
          const size_t t = std::min<size_t>(l, max_len);
          ret.write<bnh>(enc.bits.code[z_code]);
          ret.write<bnh>({b, t - min_len});
          l -= t;
        }
        return l;
      };

      ret.write<bnh>({9, table.size()});
      for (size_t i = 0; i < table.size(); ) {
        if (table[i] == 0) {
          size_t l = encode::run_length(table, i, 0);
          i += l;
          l = write_zeros(l, 2, 9, 0x14);
          l = write_zeros(l, 1, 4, 0x03);
          l = write_zeros(l, 0, 0, 0x01);
          assert(l == 0);
        } else {
          if (table[i] > config.max_bitlen[h_tag::code]) {
            // this should not happen.
            throw std::logic_error("This algorithm cannot compress the given data.");
          }
          ret.write<bnh>(enc.bits.code[table[i] + 2]);
          i += 1;
        }
      }
    };

    ret.write<bnh>({16, best_commands.size()});
    write_huff_bits(best_enc.bits, 5, 3, 0x13, h_tag::bits);
    write_table(best_enc);
    write_huff_bits(best_enc.ofs, 4, -1, 0x0e, h_tag::lz_ofs);

    size_t adr = 0;
    for (const auto& cmd : best_commands) {
      const auto [tag, oi] = cmd.type;
      switch (tag) {
      case uncomp: {
        ret.write<bnh>(best_enc.code.code[input[begin + adr]]);
      } break;
      case lz: {
        ret.write<bnh>(best_enc.code.code[cmd.len + huff_offset]);
        const auto& o = lz_ofs[oi];
        const size_t d = adr - cmd.lz_ofs();
        assert(o.min <= d && d <= o.max);
        ret.write<bnh>(best_enc.ofs.code[oi]);
        ret.write<bnh>({o.bitlen, d - o.min});
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    assert(adr == size);
    begin = end;
  }
  return ret.out;
}

} // namespace

std::vector<uint8_t> doraemon_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  const std::string filename = "a.bin"; // dummy filename
  const size_t header_size = 0x16 + filename.size() + 2 + 1;

  const lha_config config {
    .header_size = header_size + 2,
    .command_limit = 0x8000,
    .use_old_encoding = true,
    .allow_empty = true,
    .lz_ofs_bits = 0x0d,
    .max_bitlen = {0x10, 0x10, 0x10},
    .allow_zero_bits = [](int tag, size_t word) { return tag != h_tag::code || (word == 0); }
  };
  const std::string comp_method = "-lh5-";

  auto ret = doraemon_comp_core(input, config);
  ret[0x00] = header_size;
  std::ranges::copy(comp_method, ret.begin() + 2);
  write32(ret, 0x07, ret.size() - (header_size + 2));
  write32(ret, 0x0b, input.size());
  ret[0x13] = 0x20;
  ret[0x14] = 0x01;
  ret[0x15] = filename.size();
  std::ranges::copy(filename, ret.begin() + 0x16);
  write16(ret, 0x16 + filename.size(), utility::crc16(input));
  ret[header_size - 1] = 0x4d;
  write16(ret, header_size, 0x0000);

  if (input.size() > 0 && input.size() <= ret.size() - (2 + header_size)) {
    ret.resize(2 + header_size + input.size());
    ret[0x05] = '0';
    write32(ret, 0x07, ret.size() - (header_size + 2));
    std::ranges::copy(input, ret.begin() + 2 + header_size);
  }

  size_t check_sum = 0;
  for (size_t i = 2; i < header_size + 2; ++i) check_sum += ret[i];
  ret[0x01] = check_sum;

  return ret;
}

std::vector<uint8_t> olivias_mystery_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);
  const lha_config config {
    .header_size = 0,
    .command_limit = 0x8000,
    .use_old_encoding = true,
    .allow_empty = true,
    .lz_ofs_bits = 0x0d,
    .max_bitlen = {0x10, 0x10, 0x10},
    .allow_zero_bits = [](int, size_t) { return true; }
  };
  return doraemon_comp_core(input, config);
}

std::vector<uint8_t> shima_kousaku_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  // Caution: The decompressor seems have a bug.
  //          It may not work if some codeword has more than 12 bits.
  const lha_config config {
    .header_size = 2,
    .command_limit = 0x8000,
    .use_old_encoding = true,
    .allow_empty = true,
    .lz_ofs_bits = 0x0e,
    .max_bitlen = {0x0c, 0x10, 0x10},
    .allow_zero_bits = [](int, size_t) { return true; }
  };
  auto ret = doraemon_comp_core(input, config);
  write16(ret, 0, input.size());
  return ret;
}

std::vector<uint8_t> yatterman_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  // Caution: The decompressor seems have a bug.
  // - It may not work if some codeword has more than 13 bits. (cf. $80:B407, $80:B40A)
  // - It may not work if some (lz-offset) codeword has more than 8 bits. (cf. $80:B4BE)
  const lha_config config = {
    .header_size = 2,
    .command_limit = 0x8000,
    .use_old_encoding = false,
    .allow_empty = true,
    .lz_ofs_bits = 0x0e,
    .max_bitlen = {0x0d, 0x08, 0x0f},
    .allow_zero_bits = [](int, size_t) { return false; }
  };
  auto ret = doraemon_comp_core(input, config);
  write16(ret, 0, input.size());
  if (ret.size() % 2 == 1) ret.push_back(0);
  for (size_t i = 2; i < ret.size(); i += 2) std::swap(ret[i], ret[i + 1]);
  return ret;
}

std::vector<uint8_t> time_cop_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);
  const lha_config config {
    .header_size = 6,
    .command_limit = 0x7fff,
    .use_old_encoding = false,
    .allow_empty = false,
    .lz_ofs_bits = 0x0d,
    .max_bitlen = {0x10, 0x10, 0x10},
    .allow_zero_bits = [](int, size_t) { return false; }
  };
  assert(config.command_limit <= 0x7fff);
  auto ret = doraemon_comp_core(input, config);
  // The data should be aligned to an even address if it crosses banks.
  write16(ret, 0, input.size() >> 0);
  write16(ret, 2, input.size() >> 16); // unknown
  write16(ret, 4, utility::crc16(input));
  return ret;
}

} // namespace sfc_comp
