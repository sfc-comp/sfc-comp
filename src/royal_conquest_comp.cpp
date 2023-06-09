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

  static constexpr auto lz_ofs_tab = to_vranges({
    {0x0001,  9, 0x00 << 1}, // [00, 1f] _
    {0x0041, 10, 0x20 << 2}, // [20, 4f] __
    {0x0101, 11, 0x50 << 3}, // [50, 8f] ___
    {0x0301, 12, 0x90 << 4}, // [90, bf] ____
    {0x0601, 13, 0xc0 << 5}, // [c0, ef] _____
    {0x0c01, 14, 0xf0 << 6}, // [f0, ff] ______
  }, 0x0fff); // not 0x1000 (cf. $02:E1AB)

  static constexpr size_t word_max_count = 0x13a;
  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = lz_min_len + 0xff;
  static constexpr size_t lz_huff_offset = 0x0100 - lz_min_len;

  static constexpr size_t pad = lz_max_len;

  std::vector<uint8_t> input(in.size() + pad, 0);
  std::ranges::copy(in, input.begin() + pad);

  const auto lz_memo = [&] {
    std::vector<std::array<encode::lz_data, lz_ofs_tab.size()>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, lz_ofs_tab, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find(i, max_ofs, lz_min_len);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  static constexpr auto lz_lens = create_array<size_t, (lz_max_len - lz_min_len) + 1>([&](size_t i) {
    return i + lz_min_len;
  });

  auto huff_bitlens = [&]{
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto h = encode::huffman(utility::freq_u8(in), true);
    for (const auto w : h.words) ret[w] = h.code[w].bitlen;
    return ret;
  }();

  const auto update_costs = [&](const encode::huffman_t& huff, size_t penalty = 0) {
    const auto& codewords = huff.code;
    size_t longest_bit_size = 0;
    for (const auto w : huff.words) longest_bit_size = std::max<size_t>(longest_bit_size, codewords[w].bitlen);
    huff_bitlens.assign(huff_bitlens.size(), longest_bit_size + penalty);
    for (const auto w : huff.words) huff_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;

  size_t penalty = ilog2(2 * input.size() + 1) + 9;
  while (true) {
    solver<tag> dp(input.size());

    for (size_t i = input.size(); i-- > pad; ) {
      dp.update_matrix(i, lz_ofs_tab, lz_lens,
        [&](size_t li) { return huff_bitlens[lz_lens[li] + lz_huff_offset]; },
        [&](size_t oi) { return lz_memo[i][oi]; },
        [&](size_t oi, size_t) -> tag { return {lz, oi}; }
      );
      dp.update(i, 1, huff_bitlens[input[i]], {uncomp, 0});
    }

    std::array<size_t, 0x100 + lz_lens.size()> code_counts = {};
    {
      size_t adr = pad;
      for (const auto& cmd : dp.optimal_path(pad)) {
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
        const auto c = huff.code[w];
        const size_t r_zero = std::min<size_t>(std::countr_zero(c.val), c.bitlen);
        tree.write<bnh>({r_zero, low_bits_mask(r_zero)});
        tree.write<b1>(false);
        ret[4 + 0x4f + (i >> 3)] |= (w >> 8) << (7 - (i & 7));
        ret[4 + 0x77 + i] = w & 0xff;
      }
      assert(tree.size() <= 0x4f);
      for (size_t i = 0; i < tree.size(); ++i) ret[4 + i] = tree[i];

      size_t adr = pad;
      for (const auto& cmd : dp.optimal_path(pad)) {
        const auto [tag, oi] = cmd.type;
        switch (tag) {
        case uncomp: {
          ret.write<bnh>(huff.code[input[adr]]);
        } break;
        case lz: {
          ret.write<bnh>(huff.code[cmd.len + lz_huff_offset]);
          const auto o = lz_ofs_tab[oi];
          ret.write<bnh>({o.bitlen, o.val + (adr - cmd.lz_ofs() - o.min)});
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
