#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

std::vector<uint8_t> royal_conquest_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);

  enum method { uncomp, lz };
  using tag = tag_o<method>;

  static constexpr auto lz_ofs_tab = std::to_array<vrange>({
    vrange(0x0001, 0x0040,  9, 0x00 << 1), // [00, 1f] _
    vrange(0x0041, 0x0100, 10, 0x20 << 2), // [20, 4f] __
    vrange(0x0101, 0x0300, 11, 0x50 << 3), // [50, 8f] ___
    vrange(0x0301, 0x0600, 12, 0x90 << 4), // [90, bf] ____
    vrange(0x0601, 0x0c00, 13, 0xc0 << 5), // [c0, ef] _____
    vrange(0x0c01, 0x0fff, 14, 0xf0 << 6), // [f0, ff] ______ // not 0x1000 (cf. $02:E1AB)
  });

  static constexpr size_t word_max_count = 0x13a;
  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = lz_min_len + 0xff;
  static constexpr size_t lz_huff_offset = 0x0100 - lz_min_len;

  static constexpr size_t pad = lz_max_len;

  std::vector<uint8_t> input(in.size() + pad, 0);
  std::copy(in.begin(), in.end(), input.begin() + pad);

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, lz_ofs_tab.size()>> lz_memo(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    size_t o = lz_ofs_tab.size() - 1;
    auto res_lz = lz_helper.find_best_closest(i, lz_ofs_tab[o].max, lz_max_len);
    if (res_lz.len < lz_min_len) res_lz = {0, 0};
    lz_memo[i][o] = res_lz;
    for (; o-- > 0; ) {
      const size_t d = i - res_lz.ofs;
      const size_t max_ofs = lz_ofs_tab[o].max;
      if (res_lz.len >= lz_min_len && d > max_ofs) {
        res_lz = lz_helper.find_best_closest(i, max_ofs, lz_max_len);
      }
      lz_memo[i][o] = res_lz;
    }
    lz_helper.add_element(i);
  }

  static constexpr auto lz_lens = create_array<size_t, (lz_max_len - lz_min_len) + 1>([&](size_t i) {
    return i + lz_min_len;
  });

  auto huff_bitlens = [&]{
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto h = encode::huffman(utility::freq_u8(in), true);
    for (const auto w : h.words) ret[w] = h.codewords[w].bitlen;
    return ret;
  }();

  const auto update_costs = [&](const encode::huffman_t& huff, size_t penalty = 0) {
    const auto& codewords = huff.codewords;
    size_t longest_bit_size = 0;
    for (const auto w : huff.words) longest_bit_size = std::max<size_t>(longest_bit_size, codewords[w].bitlen);
    huff_bitlens.assign(huff_bitlens.size(), std::max(penalty, longest_bit_size + 9));
    for (const auto w : huff.words) huff_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;

  size_t penalty = ilog2(2 * input.size() + 1) + 9;
  while (true) {
    sssp_solver<tag> dp(input.size(), pad);

    for (size_t i = pad; i < input.size(); ++i) {
      dp.update(i, 1, 1, [&](size_t) { return huff_bitlens[input[i]]; }, {uncomp, 0});
      for (size_t o = lz_ofs_tab.size(); o-- > 0; ) {
        const auto& res_lz = lz_memo[i][o];
        if (res_lz.len < lz_min_len) break;
        if ((i - res_lz.ofs) < lz_ofs_tab[o].min) continue;
        dp.update_lz_table(i, lz_lens, res_lz, [&](size_t j) {
          return huff_bitlens[lz_lens[j] + lz_huff_offset] + lz_ofs_tab[o].bitlen;
        }, {lz, o});
      }
    }

    const auto commands = dp.commands(pad);

    std::array<size_t, 0x100 + lz_lens.size()> code_counts = {};
    {
      size_t adr = pad;
      for (const auto& cmd : commands) {
        switch (cmd.type.tag) {
        case uncomp: code_counts[input[adr]] += 1; break;
        case lz: code_counts[cmd.len + lz_huff_offset] += 1; break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());
    }

    const auto huff = encode::huffman(code_counts, true);

    if (huff.words.size() > word_max_count) {
      size_t excess = huff.words.size() - word_max_count;
      for (ptrdiff_t i = huff.words.size() - 1; i >= 0; --i) {
        size_t w = huff.words[i];
        if (w < 0x100) continue;
        code_counts[w] = 0;
        --excess;
        if (excess == 0) break;
      }
      const auto huff_fixed = encode::huffman(code_counts, true);
      update_costs(huff_fixed, penalty);
      penalty = penalty * 1.25 + 1;
      if (penalty >= 1e9) {
        // this should not happen.
        throw std::runtime_error("Failed to compress the given input.");
      }
    } else {
      if (huff.words.size() == 0) {
        throw std::logic_error("Input is empty.");
      }
      using namespace data_type;

      writer_b8_l ret(4 + 0x77 + word_max_count);
      writer_b8_h tree;
      for (size_t i = 0; i < huff.words.size(); ++i) {
        const auto w = huff.words[i];
        const auto c = huff.codewords[w];
        const size_t r_zero = std::min<size_t>(std::countr_zero(c.val), c.bitlen);
        tree.write<bnh>({r_zero, (size_t(1) << r_zero) - 1});
        tree.write<b1>(false);
        ret[4 + 0x4f + (i >> 3)] |= (w >> 8) << (7 - (i & 7));
        ret[4 + 0x77 + i] = w & 0xff;
      }
      assert(tree.size() <= 0x4f);
      for (size_t i = 0; i < tree.size(); ++i) ret[4 + i] = tree[i];

      size_t adr = pad;
      for (const auto& cmd : commands) {
        switch (cmd.type.tag) {
        case uncomp: {
          ret.write<bnh>(huff.codewords[input[adr]]);
        } break;
        case lz: {
          ret.write<bnh>(huff.codewords[cmd.len + lz_huff_offset]);
          const auto o = lz_ofs_tab[cmd.type.oi];
          ret.write<bnh>({o.bitlen, o.val + (adr - cmd.lz_ofs - o.min)});
        } break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      write32(ret.out, 0, in.size());
      assert(adr == input.size());

      if (best.empty() || ret.size() < best.size()) {
        best = std::move(ret.out);
      } else {
        break;
      }
      update_costs(huff);
    }
  }

  return best;
}

} // namespace sfc_comp
