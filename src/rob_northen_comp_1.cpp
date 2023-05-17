#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

std::vector<uint8_t> rob_northen_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x100000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr size_t lz_min_len = 2;
  static constexpr size_t lz_ofs_max_bits = 15;
  static constexpr size_t lz_len_max_bits = 15; // <= 15
  static constexpr size_t ulen_max_bits = 15; // <= 15

  struct rnc1_huff {
    encode::huffman_t ulen;
    encode::huffman_t ofs;
    encode::huffman_t len;
  };

  static constexpr auto lz_offsets = create_array<vrange, lz_ofs_max_bits + 1>([](size_t k) {
    if (k == 0) return vrange(1, 1, 0, 0);
    return vrange((1 << (k - 1)) + 1, 1 << k, k - 1, 0);
  });
  static constexpr auto lz_lens = create_array<vrange, lz_len_max_bits + 1>([](size_t k) {
    if (k == 0) return vrange(lz_min_len, lz_min_len, 0, 0);
    return vrange((1 << (k - 1)) + lz_min_len, (1 << k) + (lz_min_len - 1), k - 1, 0);
  });
  static constexpr auto ulens = create_array<vrange, ulen_max_bits + 1>([](size_t k) {
    if (k == 0) return vrange(0, 0, 0, 0);
    return vrange(1 << (k - 1), (1 << k) - 1, k - 1, 0);
  });
  static constexpr size_t chunk_size = ulens.back().max;

  const auto max_elem = [&](std::span<const size_t> words) -> ptrdiff_t {
    if (words.size() == 0) return -1; // This value should be -1.
    return *std::max_element(words.begin(), words.end());
  };

  const auto lz_memo = [&] {
    std::vector<std::array<encode::lz_data, lz_offsets.size()>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, lz_offsets, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find(i, max_ofs, lz_min_len);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  using namespace data_type;
  writer_b16_l ret(0x12); ret.write<bnl>({2, 0});

  size_t chunk = 0;
  for (size_t begin = 0; begin < input.size(); ) {
    // [Todo] Find a reasonable initialization.
    std::vector<vrange> curr_uncomp(ulens.begin(), ulens.end());
    for (size_t k = 0; k < ulens.size(); ++k) curr_uncomp[k].bitlen += k + 1;

    std::vector<vrange> curr_lz_len(lz_lens.begin(), lz_lens.end());
    for (size_t k = 0; k < lz_lens.size(); ++k) curr_lz_len[k].bitlen += k + 1;

    std::vector<vrange> curr_lz_ofs(lz_offsets.begin(), lz_offsets.end());
    for (size_t k = 0; k < lz_offsets.size(); ++k) curr_lz_ofs[k].bitlen += 4;

    const size_t end = std::min(input.size(), begin + chunk_size);
    const size_t size = end - begin;

    using node_type = solver<tag>::node;
    size_t best_cost = std::numeric_limits<size_t>::max();
    std::vector<node_type> best_commands;
    rnc1_huff best_huff;

    for (bool updated = true; updated; ) {
      const auto shift_lz = [&](const encode::lz_data& p) -> encode::lz_data {
        return {p.ofs - begin, p.len}; // Note: can be negative
      };

      solver<tag> dp0(size, -1), dp1(size);
      auto c0_0 = dp0.c<0>(lz_lens.back().max);
      auto c8_1 = dp1.c<8>(ulens.back().max);

      for (size_t i = size; ; ) {
        dp0.update_c(i, 0, dp1[i].cost + curr_uncomp[0].bitlen, {uncomp, 0, 0});
        c0_0.update(i);
        if (i-- == 0) break;
        dp0.update(i, std::span(curr_uncomp.begin() + 1, curr_uncomp.end()), c8_1, 0,
          [&](size_t li) -> tag { return {uncomp, 0, li + 1}; });
        dp1.update_matrix(i, curr_lz_ofs, curr_lz_len, c0_0, 0,
          [&](size_t oi) { return shift_lz(lz_memo[begin + i][oi]); },
          [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
        );
        c8_1.update(i);
      }

      if (dp0.optimal_cost() == dp0.infinite_cost) {
        throw std::logic_error("Failed to compress. chunk_size may be too large.");
      }

      struct counter_t {
        std::array<size_t, ulen_max_bits + 1> ulen = {};
        std::array<size_t, lz_len_max_bits + 1> len = {};
        std::array<size_t, lz_ofs_max_bits + 1> ofs = {};
      } counter;

      const auto commands = [&] {
        std::vector<node_type> ret;
        size_t adr = 0;
        for (size_t ty = 0; adr < size || (adr == size && ty == 0); ty ^= 1) {
          auto cmd = (ty == 0) ? dp0[adr] : dp1[adr];
          const auto [tag, oi, li] = cmd.type;
          if (tag == uncomp) counter.ulen[li] += 1;
          else counter.ofs[oi] += 1, counter.len[li] += 1;
          adr += cmd.len;
          ret.emplace_back(cmd);
        }
        assert(adr == size);
        return ret;
      }();

      const auto enc_huffman = [&](std::span<const size_t> counts) {
        auto ret = encode::huffman(counts, true);
        // 0-bit codes are not allowed.
        if (ret.words.size() == 1) ret.code[ret.words[0]].bitlen += 1;
        return ret;
      };

      const rnc1_huff huff(
        enc_huffman(counter.ulen),
        enc_huffman(counter.ofs),
        enc_huffman(counter.len)
      );

      const auto calc_cost = [&](const rnc1_huff& huff,
          std::span<const node_type> commands, const counter_t& counter) -> size_t {
        size_t cost = 0;
        cost += 5 + 4 * (1 + max_elem(huff.ulen.words));
        cost += 5 + 4 * (1 + max_elem(huff.ofs.words));
        cost += 5 + 4 * (1 + max_elem(huff.len.words));
        cost += 16;
        for (const auto& cmd : commands) if (cmd.type.tag == uncomp) cost += cmd.len * 8;
        for (size_t w : huff.ulen.words) cost += counter.ulen[w] * (huff.ulen.code[w].bitlen + ulens[w].bitlen);
        for (size_t w : huff.ofs.words) cost += counter.ofs[w] * (huff.ofs.code[w].bitlen + lz_offsets[w].bitlen);
        for (size_t w : huff.len.words) cost += counter.len[w] * (huff.len.code[w].bitlen + lz_lens[w].bitlen);
        return cost;
      };

      const auto update_costs = [&](const rnc1_huff& huff) {
        const auto update = [&](const encode::huffman_t& huff, std::span<const vrange> base,
            std::vector<vrange>& costs, size_t penalty = 1) {
          size_t longest_bitlen = 0;
          for (const auto w : huff.words) longest_bitlen = std::max<size_t>(longest_bitlen, huff.code[w].bitlen);
          for (auto& v : costs) v.bitlen = longest_bitlen + penalty;
          for (const auto w : huff.words) costs[w].bitlen = huff.code[w].bitlen + base[w].bitlen;
        };
        update(huff.ulen, ulens, curr_uncomp);
        update(huff.ofs, lz_offsets, curr_lz_ofs);
        update(huff.len, lz_lens, curr_lz_len);
      };

      const size_t chunk_cost = calc_cost(huff, commands, counter);
      if (chunk_cost < best_cost) {
        best_cost = chunk_cost;
        best_commands = std::move(commands);
        best_huff = std::move(huff);
        update_costs(best_huff);
      } else {
        updated = false;
      }
    }

    const auto write_bitlens = [&](const encode::huffman_t& huff) {
      const size_t total_count = max_elem(huff.words) + 1;
      ret.write<bnl>({5, total_count});
      for (size_t i = 0; i < total_count; ++i) {
        ptrdiff_t bits = huff.code[i].bitlen;
        if (bits <= 0) bits = 0;
        // this would not happen.
        if (bits >= 16) throw std::runtime_error("Failed to compress the given data. (bits >= 16)");
        ret.write<bnl>({4, size_t(bits)});
      }
    };

    const auto prev_cost = ret.bit_length();
    write_bitlens(best_huff.ulen);
    write_bitlens(best_huff.ofs);
    write_bitlens(best_huff.len);
    assert(best_commands.size() % 2 == 1);
    const size_t command_count = (best_commands.size() + 1) / 2;
    if (command_count >= 0x10000) throw std::logic_error("chunk_size may be too large.");
    ret.write<bnl>({16, command_count});

    size_t adr = 0;
    for (const auto& cmd : best_commands) {
      const auto [tag, oi, li] = cmd.type;
      switch (tag) {
      case uncomp: {
        ret.write<bnh>(best_huff.ulen.code[li]);
        if (cmd.len > 0) {
          assert(ulens[li].min <= cmd.len && cmd.len <= ulens[li].max);
          ret.write<bnl>({ulens[li].bitlen, cmd.len - ulens[li].min});
          ret.write<d8n>({cmd.len, &input[begin + adr]});
        }
      } break;
      case lz: {
        const auto& o = lz_offsets[oi];
        const size_t d = (adr - cmd.lz_ofs());
        assert(o.min <= d && d <= o.max);
        ret.write<bnh>(best_huff.ofs.code[oi]);
        ret.write<bnl>({o.bitlen, d - o.min});

        const auto l = lz_lens[li];
        assert(l.min <= cmd.len && cmd.len <= l.max);
        ret.write<bnh>(best_huff.len.code[li]);
        ret.write<bnl>({l.bitlen, cmd.len - l.min});
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    assert(best_cost == ret.bit_length() - prev_cost);
    assert(adr == size);
    begin = end;

    ++chunk;
    if (chunk == 0x100) throw std::runtime_error("This algorithm cannot compress the given data. (chunk >= 0x100)");
  }

  write32b(ret.out, 0, 0x524e4301);
  write32b(ret.out, 4, input.size());
  write32b(ret.out, 8, ret.size() - 18);
  write16b(ret.out, 12, utility::crc16(input, 0, input.size()));
  write16b(ret.out, 14, utility::crc16(ret.out, 18, ret.size() - 18));
  ret[0x10] = 0; // leeway
  ret[0x11] = chunk;

  return ret.out;
}

} // namespace sfc_comp
