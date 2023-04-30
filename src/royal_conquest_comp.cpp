#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

std::vector<uint8_t> royal_conquest_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);

  enum Tag {
    uncomp, lz
  };

  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return true;
      return true;
    }
    Tag tag;
    size_t ofs_no;
  };

  struct tup {
    size_t min;
    size_t val;
    size_t bits;
  };

  static constexpr size_t lz_ofs_total = 6;
  static constexpr std::array<tup, lz_ofs_total + 1> lz_ofs_table = {
    tup({0x0001, 0x00 << 1, 9}),  // [00, 1f] _
    tup({0x0041, 0x20 << 2, 10}), // [20, 4f] __
    tup({0x0101, 0x50 << 3, 11}), // [50, 8f] ___
    tup({0x0301, 0x90 << 4, 12}), // [90, bf] ____
    tup({0x0601, 0xc0 << 5, 13}), // [c0, ef] _____
    tup({0x0c01, 0xf0 << 6, 14}), // [f0, ff] ______
    tup({0x1000, 0, 0}) // not 0x1001 (cf. $02:E1AB)
  };

  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = lz_min_len + 0xff;
  static constexpr size_t lz_huff_offset = 0x0100 - lz_min_len;
  static constexpr size_t word_max_count = 0x13a;

  static constexpr size_t pad = lz_max_len;

  std::vector<uint8_t> input(in.size() + pad, 0);
  std::copy(in.begin(), in.end(), input.begin() + pad);

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, lz_ofs_total>> lz_memo(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    size_t o = lz_ofs_total - 1;
    auto res_lz = lz_helper.find_best_closest(i, lz_ofs_table[o + 1].min - 1, lz_max_len);
    if (res_lz.len < lz_min_len) res_lz = {0, 0};
    lz_memo[i][o] = res_lz;
    for (; o-- > 0; ) {
      const size_t d = i - res_lz.ofs;
      const size_t max_ofs = lz_ofs_table[o + 1].min - 1;
      if (res_lz.len >= lz_min_len && d > max_ofs) {
        res_lz = lz_helper.find_best_closest(i, max_ofs, lz_max_len);
      }
      lz_memo[i][o] = res_lz;
    }
    lz_helper.add_element(i);
  }

  std::vector<size_t> lz_lens(0x100);
  std::iota(lz_lens.begin(), lz_lens.end(), lz_min_len);

  auto huff_bitlens = [&]{
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto codewords = encode::huffman(utility::freq_u8(in), true).codewords;
    for (size_t i = 0; i < codewords.size(); ++i) {
      if (codewords[i].bitlen >= 0) ret[i] = codewords[i].bitlen;
    }
    return ret;
  }();

  auto update_huffman_costs = [&] (const encode::huffman_result& huff, size_t penalty = 0) {
    const auto& codewords = huff.codewords;
    size_t longest_bit_size = 0;
    for (const auto w : huff.words) longest_bit_size = std::max<size_t>(longest_bit_size, codewords[w].bitlen);
    huff_bitlens.assign(huff_bitlens.size(), std::max(penalty, longest_bit_size + 9));
    for (const auto w : huff.words) huff_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;

  size_t penalty = ilog2(2 * input.size() + 1) + 9;
  while (true) {
    sssp_solver<CompType> dp(input.size(), pad);

    for (size_t i = pad; i < input.size(); ++i) {
      dp.update(i, 1, 1, [&](size_t) { return huff_bitlens[input[i]]; }, {uncomp, 0});
      for (ptrdiff_t o = lz_ofs_total - 1; o >= 0; --o) {
        const auto& res_lz = lz_memo[i][o];
        if (res_lz.len < lz_min_len) break;
        if ((i - res_lz.ofs) < lz_ofs_table[o].min) continue;
        const size_t ofs_bits = lz_ofs_table[o].bits;
        dp.update_lz_table(i, lz_lens, res_lz, [&](size_t j) {
          return huff_bitlens[lz_lens[j] + lz_huff_offset] + ofs_bits;
        }, {lz, size_t(o)});
      }
    }

    const auto commands = dp.commands(pad);

    std::array<size_t, 0x200> code_counts = {};
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
      update_huffman_costs(huff_fixed, penalty);
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
          const auto tp = lz_ofs_table[cmd.type.ofs_no];
          ret.write<bnh>({tp.bits, tp.val + (adr - cmd.lz_ofs - tp.min)});
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
      update_huffman_costs(huff);
    }
  }

  return best;
}

} // namespace sfc_comp
