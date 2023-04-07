#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

struct CostType {
  explicit constexpr CostType() : CostType(0, 0) {}
  explicit constexpr CostType(size_t cost, size_t count = 0) : cost(cost), count(count) {}

  constexpr bool operator == (const CostType& rhs) const { return cost == rhs.cost && count == rhs.count; }
  constexpr bool operator != (const CostType& rhs) const { return !(*this == rhs); }
  constexpr CostType operator + (const CostType& rhs) const {
    return CostType(cost + rhs.cost, count + rhs.count);
  }
  constexpr CostType operator + (const size_t& c) const {
    return CostType(cost + c, count + 1);
  }
  constexpr bool operator >= (const CostType& rhs) const {
    return cost >= rhs.cost;
  }
  constexpr bool operator < (const CostType& rhs) const {
    return !(*this >= rhs);
  }
  size_t cost;
  size_t count;
};

} // namespace

template <>
constexpr CostType default_cost() { return CostType(default_cost<size_t>()); }

struct lha_config {
  bool use_old_encoding = true;
  size_t header_size = 0;
  size_t command_limit = 0x4000;
  size_t code_max_bits = 0x10;
  size_t lz_ofs_bits = 0x0d;
};

template <typename Func>
std::vector<uint8_t> doraemon_comp_core(std::span<const uint8_t> input, const lha_config& config, Func use_zero_bits) {
  enum Tag {
    uncomp, lz
  };

  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return true;
      return ofs_no == rhs.ofs_no;
    }
    Tag tag;
    size_t ofs_no;
  };

  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = 0x100;

  static constexpr size_t lz_huff_offset = 0x0100 - lz_min_len;
  static constexpr size_t lz_ofs_buff_bits = 0x0f;
  assert(config.lz_ofs_bits <= lz_ofs_buff_bits);

  assert(1 <= config.command_limit && config.command_limit <= 0xffff);

  std::vector<size_t> lz_ofs_table(config.lz_ofs_bits + 1);
  for (size_t i = 0; i < lz_ofs_table.size(); ++i) lz_ofs_table[i] = 1 << i;

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, lz_ofs_buff_bits + 1>> lz_memo(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    auto res_lz = lz_helper.find_best_closest(i, lz_ofs_table.back(), lz_max_len);
    if (res_lz.len < lz_min_len) res_lz = {0, 0};
    size_t o = config.lz_ofs_bits;
    lz_memo[i][o] = res_lz;
    for (; o-- > 0; ) {
      const size_t d = i - res_lz.ofs;
      if (res_lz.len >= lz_min_len && d > lz_ofs_table[o]) {
        res_lz = lz_helper.find_best_closest(i, lz_ofs_table[o], lz_max_len);
      }
      lz_memo[i][o] = res_lz;
    }
    lz_helper.add_element(i);
  }

  std::vector<size_t> lz_lens(lz_max_len - lz_min_len + 1);
  std::iota(lz_lens.begin(), lz_lens.end(), lz_min_len);

  constexpr auto bit_cost = [](size_t b) -> size_t {
    return (b == 0) ? 0 : (b - 1);
  };

  using namespace data_type;
  writer_b ret;
  for (size_t i = 0; i < config.header_size; ++i) ret.write<d8>(0);

  using command_type = sssp_solver<CompType, CostType>::vertex_type;
  sssp_solver<CompType, CostType> dp(input.size());

  size_t begin = 0;
  size_t last_index = begin;

  while (begin < input.size()) {
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> huff_bitsizes(0x100 + (lz_max_len - lz_min_len) + 1, 9);
    for (size_t i = 0x0100; i < huff_bitsizes.size(); ++i) huff_bitsizes[i] = 18;
    std::vector<size_t> huff_lz_ofs_bits(config.lz_ofs_bits + 1, 4);

    struct encodes {
      encode::huffman_result code;
      encode::huffman_result lz_ofs;
      encode::huffman_result bits;
    };

    auto update_costs = [&] (const encode::huffman_result& huff, std::vector<size_t>& costs, size_t penalty = 9) {
      size_t worst_bit_count = 0;
      for (const auto w : huff.words) worst_bit_count = std::max<size_t>(worst_bit_count, huff.codewords[w].bit_count);
      costs.assign(costs.size(), worst_bit_count + penalty);
      for (const auto w : huff.words) costs[w] = huff.codewords[w].bit_count;
    };

    auto update_huffman_costs = [&] (const encodes& enc) {
      update_costs(enc.code, huff_bitsizes);
      update_costs(enc.lz_ofs, huff_lz_ofs_bits);
    };

    size_t best_cost = std::numeric_limits<size_t>::max();
    size_t best_last_index = 0;
    std::vector<command_type> best_commands;
    encodes best_enc;

    while (true) {
      dp[begin].cost = CostType(0);
      size_t e = std::min(input.size(), last_index + lz_max_len);
      for (size_t i = begin + 1; i <= e; ++i) dp[i].cost = dp.default_cost;

      size_t index = begin;
      for (; index < input.size(); ++index) {
        if (dp[index].cost.count >= config.command_limit) break;
        dp.update(index, 1, 1, [&](size_t) { return huff_bitsizes[input[index]]; }, {uncomp, 0});
        for (ptrdiff_t o = config.lz_ofs_bits; o >= 0; --o) {
          const auto& res_lz = lz_memo[index][o];
          if (res_lz.len < lz_min_len) break;
          const size_t lz_min_dist = (o == 0) ? lz_ofs_table[0] : (lz_ofs_table[o - 1] + 1);
          if ((index - res_lz.ofs) < lz_min_dist) continue;
          const size_t ofs_bits = huff_lz_ofs_bits[o] + bit_cost(o);
          dp.update_lz_table(index, lz_lens, res_lz, [&](size_t j) {
            return huff_bitsizes[lz_lens[j] + lz_huff_offset] + ofs_bits;
          }, {lz, size_t(o)});
        }
      }
      last_index = index;

      struct Counter {
        Counter(size_t lz_ofs_bits) : code(0x100 + (lz_max_len - lz_min_len) + 1, 0),
                                      lz_ofs(lz_ofs_bits + 1, 0) {}
        std::vector<size_t> code;
        std::vector<size_t> lz_ofs;
      } counter(config.lz_ofs_bits);

      const std::vector<command_type> commands = [&] {
        std::vector<command_type> ret;
        size_t adr = index;
        while (adr > begin) {
          auto cmd = dp[adr];
          adr -= cmd.len;
          switch (cmd.type.tag) {
          case uncomp: {
            counter.code[input[adr]] += 1;
          } break;
          case lz: {
            counter.code[cmd.len + lz_huff_offset] += 1;
            counter.lz_ofs[cmd.type.ofs_no] += 1;
          } break;
          default: assert(0);
          }
          ret.emplace_back(cmd);
        }
        assert(adr == begin);
        std::reverse(ret.begin(), ret.end());
        return ret;
      }();

      const auto encode_codes = [&](std::span<const size_t> counts, bool codes) {
        auto ret = encode::huffman(counts, true);
        if (ret.words.size() == 1) {
          if (!use_zero_bits(codes, ret.words[0])) {
            ret.codewords[ret.words[0]].bit_count += 1;
          }
        }
        return ret;
      };

      const auto max_elem = [&](std::span<const size_t> words) -> ptrdiff_t {
        if (words.size() == 0) return -1; // This value must be -1.
        return *std::max_element(words.begin(), words.end());
      };

      const auto bits_table = [&](const encode::huffman_result& huff) {
        const size_t table_size = max_elem(huff.words) + 1;
        std::vector<uint8_t> bits(table_size, 0);
        for (const auto w : huff.words) bits[w] = huff.codewords[w].bit_count;
        return bits;
      };

      const auto encode_table = [&](const encode::huffman_result& huff) {
        const auto bits = bits_table(huff);
        if (bits.empty()) return encode::huffman_result();

        const size_t max_bit_size = huff.codewords[huff.words.back()].bit_count;
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
        return encode_codes(counts, false);
      };

      const auto calc_cost = [&](const encodes& enc, std::span<const command_type>, const Counter& counter) -> size_t {
        size_t ret = 0;
        for (const auto& w : enc.code.words) ret += enc.code.codewords[w].bit_count * counter.code[w];
        for (const auto& w : enc.lz_ofs.words) ret += (enc.lz_ofs.codewords[w].bit_count + bit_cost(w)) * counter.lz_ofs[w];
        return ret;
      };

      encodes enc;
      enc.code = encode_codes(counter.code, true);
      enc.lz_ofs = encode_codes(counter.lz_ofs, false);
      enc.bits = encode_table(enc.code);

      const double estimated_cost = calc_cost(enc, commands, counter);

      if (estimated_cost < best_cost) {
        best_cost = estimated_cost;
        best_commands = std::move(commands);
        best_enc = std::move(enc);
        best_last_index = last_index;
        update_huffman_costs(best_enc);
      } else {
        const auto write_huff_bits =
            [&](const encode::huffman_result& huff, const size_t s, const size_t t, const size_t lim) {
          auto bits = bits_table(huff);
          if (bits.size() > lim) std::runtime_error("This algorithm cannot compress the given data.");

          if (config.use_old_encoding) {
            if (huff.words.size() == 0) {
              ret.write<b8hn_h>({s, 1});
              ret.write<b8hn_h>({3, 0});
              return;
            } else if (huff.words.size() == 1) {
              // Caution: These codes might not work on the original LHA algorithm.
              ret.write<b8hn_h>({s, 0});
              ret.write<b8hn_h>({s, huff.words[0]});
              return;
            }
          }

          ret.write<b8hn_h>({s, bits.size()});
          for (size_t i = 0; i < bits.size(); ) {
            if (bits[i] < 7) {
              ret.write<b8hn_h>({3, bits[i]});
            } else {
              ret.write<b8hn_h>({3, 7});
              for (size_t j = bits[i]; j >= 8; --j) ret.write<b8hn_h>({1, 1});
              ret.write<b8hn_h>({1, 0});
            }
            i += 1;
            if (i == t) {
              size_t l = encode::run_length(bits, i, 0);
              l = std::min<size_t>(l, 3);
              if (bits[i] != 0) l = 0;
              ret.write<b8hn_h>({2, l});
              i += l;
            }
          }
        };

        const auto write_table = [&](const encodes& enc) {
          const auto table = bits_table(enc.code);

          const auto& words = enc.code.words;

          if (config.use_old_encoding) {
            if (words.size() == 0) {
              std::logic_error("This should not happen");
            }
            if (words.size() == 1) {
              if (use_zero_bits(true, words[0])) {
                // Caution: These codes might not work on the original LHA algorithm.
                ret.write<b8hn_h>({9, 0});
                ret.write<b8hn_h>({9, words[0]});
                return;
              }
            }
          }

          ret.write<b8hn_h>({9, table.size()});

          const auto write_zeros = [&](size_t l, const size_t z_code, const size_t b, const size_t min_len) -> size_t {
            const size_t max_len = min_len + (1 << b) - 1;
            while (l >= min_len) {
              const size_t t = std::min<size_t>(l, max_len);
              const auto c = enc.bits.codewords[z_code];
              ret.write<b8hn_h>({size_t(c.bit_count), c.val});
              if (b > 0) ret.write<b8hn_h>({b, t - min_len});
              l -= t;
            }
            return l;
          };

          for (size_t i = 0; i < table.size(); ) {
            if (table[i] == 0) {
              size_t l = encode::run_length(table, i, 0);
              i += l;
              l = write_zeros(l, 2, 9, 0x14);
              l = write_zeros(l, 1, 4, 0x03);
              l = write_zeros(l, 0, 0, 0x01);
              assert(l == 0);
            } else {
              if (table[i] > config.code_max_bits) {
                throw std::runtime_error("This algorithm cannot compress the given data.");
              }
              const auto c = enc.bits.codewords[table[i] + 2];
              ret.write<b8hn_h>({size_t(c.bit_count), c.val});
              i += 1;
            }
          }
        };

        ret.write<b8hn_h>({16, commands.size()});
        write_huff_bits(best_enc.bits, 5, 3, 0x13);
        write_table(best_enc);
        write_huff_bits(best_enc.lz_ofs, 4, -1, 0x0e);

        size_t adr = begin;
        for (const auto& cmd : best_commands) {
          switch (cmd.type.tag) {
          case uncomp: {
            const auto c = best_enc.code.codewords[input[adr]];
            ret.write<b8hn_h>({size_t(c.bit_count), c.val});
          } break;
          case lz: {
            const auto cl = best_enc.code.codewords[cmd.len + lz_huff_offset];
            ret.write<b8hn_h>({size_t(cl.bit_count), cl.val});

            const size_t o = cmd.type.ofs_no;
            const size_t d = adr - cmd.lz_ofs;
            const size_t min_ofs = (o == 0) ? 1 : lz_ofs_table[o - 1] + 1;
            assert(min_ofs <= d && d <= lz_ofs_table[o]);

            const auto co = best_enc.lz_ofs.codewords[o];
            ret.write<b8hn_h>({size_t(co.bit_count), co.val});
            if (o >= 2) ret.write<b8hn_h>({o - 1, d - min_ofs});
          } break;
          default: assert(0);
          }
          adr += cmd.len;
        }
        assert(adr == best_last_index);
        begin = adr;
        break;
      }
    }
  }

  return ret.out;
}

std::vector<uint8_t> doraemon_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  const std::string filename = "a.bin"; // dummy filename
  const size_t header_size = 0x16 + filename.size() + 2 + 1;

  lha_config config;
  config.header_size = header_size + 2;
  config.code_max_bits = 0x10;
  config.lz_ofs_bits = 0x0d;
  config.command_limit = 0x8000;
  config.use_old_encoding = true;

  auto ret = doraemon_comp_core(input, config, [](bool codes, size_t word) {
    return !codes || (word == 0);
  });

  const std::string comp_method = "-lh5-";

  ret[0x00] = header_size;
  std::copy(comp_method.begin(), comp_method.end(), ret.begin() + 2);
  write32(ret, 0x07, ret.size() - (header_size + 2));
  write32(ret, 0x0b, input.size());
  ret[0x13] = 0x20;
  ret[0x14] = 0x01;
  ret[0x15] = filename.size();
  std::copy(filename.begin(), filename.end(), ret.begin() + 0x16);
  write16(ret, 0x16 + filename.size(), utility::crc16(input));
  ret[header_size - 1] = 0x4d;
  write16(ret, header_size, 0x0000);

  if (input.size() > 0 && input.size() <= ret.size() - (2 + header_size)) {
    ret.resize(2 + header_size + input.size());
    ret[0x05] = '0';
    write32(ret, 0x07, ret.size() - (header_size + 2));
    std::copy(input.begin(), input.end(), ret.begin() + 2 + header_size);
  }

  size_t check_sum = 0;
  for (size_t i = 2; i < header_size + 2; ++i) check_sum += ret[i];
  ret[0x01] = check_sum;

  return ret;
}

std::vector<uint8_t> shima_kousaku_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  lha_config config;
  config.header_size = 2;
  config.use_old_encoding = true;

  config.lz_ofs_bits = 0x0e;

  // Caution: The decompressor seems have a bug.
  //          It may not work if some codeword has more than 12 bits.
  config.code_max_bits = 0x0c;

  // [TODO]
  // This limit is not enough to avoid runtime errors.
  config.command_limit = 0x1000;

  auto ret = doraemon_comp_core(input, config, [](bool, size_t) {
    return true;
  });
  write16(ret, 0, input.size());
  return ret;
}

std::vector<uint8_t> yatterman_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  lha_config config;
  config.header_size = 2;
  config.use_old_encoding = false;

  config.lz_ofs_bits = 0x0e;

  // Caution: The decompressor seems have a bug.
  //          It may not work if some codeword has more than 13 bits. (cf. $80:b407, $80:b40a)
  config.code_max_bits = 0x0d;

  // [TODO]
  // This limit is not enough to avoid runtime errors.
  config.command_limit = 0x2000;

  auto ret = doraemon_comp_core(input, config, [](bool, size_t) {
    return false;
  });
  write16(ret, 0, input.size());
  if (ret.size() % 2 == 1) ret.push_back(0);
  for (size_t i = 2; i < ret.size(); i += 2) std::swap(ret[i], ret[i + 1]);
  return ret;
}

} // namespace sfc_comp
