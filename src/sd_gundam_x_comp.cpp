#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

struct huffman {
  std::vector<size_t> words;
  std::vector<encode::codeword> codewords;
  std::vector<encode::node> table;
};

huffman huffman_encode(std::span<const size_t> counts) {
  auto huff = encode::huffman(counts, false);
  if (huff.words.empty()) return huffman();

  size_t depth = huff.codewords[huff.words.back()].bitlen;

  std::array<size_t, 64> depth_counts = {};
  for (const auto word : huff.words) depth_counts[huff.codewords[word].bitlen] += 1;
  std::reverse(huff.words.begin(), huff.words.end());

  huffman ret;
  std::vector<encode::node> huff_table;

  if (depth > 0) {
    size_t prev_count = 0, node_id = 0, word_id = 0;
    for (; depth > 0; --depth) {
      size_t cnt = depth_counts[depth];
      size_t next_node_id = huff_table.size();
      node_id += huff.words.size();
      for (size_t i = 0; i < (prev_count / 2); ++i) {
        huff_table.push_back({node_id, node_id + 1});
        node_id += 2;
      }
      if (prev_count & 1) {
        assert(cnt & 1);
        huff_table.push_back({node_id, word_id});
        node_id += 1; word_id += 1;
      }
      for (size_t i = 0; i < (cnt / 2); ++i) {
        huff_table.push_back({word_id, word_id + 1});
        word_id += 2;
      }
      prev_count = (cnt + prev_count) / 2;
      node_id = next_node_id;
    }
    assert(prev_count == 1);
  }

  ret.codewords = std::move(huff.codewords);
  ret.words = std::move(huff.words);
  ret.table = std::move(huff_table);
  return ret;
}

} // namespace

std::vector<uint8_t> sd_gundam_x_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  if (input.size() == 0) {
    return std::vector<uint8_t>(2, 0);
  }

  enum CompType {
    uncomp, lzs, lzl
  };

  lz_helper lz_helper(input);
  std::vector<encode::lz_data> lz_memo_s(input.size());
  std::vector<encode::lz_data> lz_memo_l(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    lz_memo_s[i] = lz_helper.find_best(i, 0x0080);
    lz_memo_l[i] = lz_helper.find_best(i, 0x0400);
    lz_helper.add_element(i);
  }

  static constexpr size_t lz_min_len = 3;
  static constexpr size_t max_code = 0x1ff;
  static constexpr size_t lz_huff_offset = 0x0100 - lz_min_len;

  std::vector<size_t> lz_lens(max_code - 0x100 + 1);
  std::iota(lz_lens.begin(), lz_lens.end(), lz_min_len);

  auto huff_bitlens = [&]{
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto codewords = encode::huffman(utility::freq_u8(input), true).codewords;
    for (size_t i = 0; i < codewords.size(); ++i) {
      if (codewords[i].bitlen >= 0) ret[i] = codewords[i].bitlen;
    }
    return ret;
  }();

  auto update_huffman_costs = [&] (const huffman& huff, size_t penalty = 0) {
    const auto& codewords = huff.codewords;
    size_t longest_bit_size = 0;
    for (const auto w : huff.words) longest_bit_size = std::max<size_t>(longest_bit_size, codewords[w].bitlen);
    huff_bitlens.assign(huff_bitlens.size(), std::max(penalty, longest_bit_size + 9));
    for (const auto w : huff.words) huff_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;

  size_t penalty = ilog2(2 * input.size() + 1) + 9;

  while (true) {
    sssp_solver<CompType> dp(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, [&](size_t) { return huff_bitlens[input[i]]; }, uncomp);
      dp.update_lz_table(i, lz_lens, lz_memo_s[i],
        [&](size_t j) { return 8 + huff_bitlens[lz_lens[j] + lz_huff_offset]; }, lzs);
      dp.update_lz_table(i, lz_lens, lz_memo_l[i],
        [&](size_t j) { return 11 + huff_bitlens[lz_lens[j] + lz_huff_offset]; }, lzl);
    }

    const auto commands = dp.commands();

    std::array<size_t, max_code + 1> code_counts = {};
    {
      size_t adr = 0;
      for (const auto& cmd : commands) {
        switch (cmd.type) {
        case uncomp: code_counts[input[adr]] += 1; break;
        case lzs:
        case lzl: code_counts[cmd.len + lz_huff_offset] += 1; break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());
      size_t non_zeros = 0;
      for (size_t i = 0; i < code_counts.size(); ++i) {
        if (code_counts[i] > 0) non_zeros += 1;
      }
      if (non_zeros == 1) code_counts[0x100] += 1;
    }

    const auto huff = huffman_encode(code_counts);

    if (huff.words.size() > 0x0180) {
      size_t excess = huff.words.size() - 0x180;
      for (size_t i = 0; i < huff.words.size(); ++i) {
        size_t w = huff.words[i];
        if (w < 0x100) continue;
        code_counts[w] = 0;
        --excess;
        if (excess == 0) break;
      }
      const auto huff_fixed = huffman_encode(code_counts);
      update_huffman_costs(huff_fixed, penalty);
      penalty = penalty * 1.25 + 1;
      if (penalty >= 1e9) {
        // this should not happen.
        throw std::runtime_error("Failed to compress the given input.");
      }
    } else {
      using namespace data_type;
      writer_b8_h ret; ret.write<bnh>({16, input.size()});

      assert(huff.words.size() >= 2);
      ret.write<bnh>({9, huff.words.size()});
      for (size_t i = 0; i < huff.words.size(); ++i) {
        ret.write<bnh>({9, huff.words[i]});
      }
      size_t prev = 0;
      auto write_child = [&](size_t v) {
        if (v == prev) {
          ret.write<b1>(false);
        } else {
          ret.write<b1, bnh>(true, {10, v});
          prev = v;
        }
        ++prev;
      };
      for (size_t i = 0; i < huff.table.size(); ++i) {
        const auto& h = huff.table[i];
        write_child(h.left);
        write_child(h.right);
      }

      size_t adr = 0;
      for (const auto& cmd : commands) {
        switch (cmd.type) {
        case uncomp: {
          const auto& c = huff.codewords[input[adr]];
          ret.write<bnh>({size_t(c.bitlen), c.val});
        } break;
        case lzs:
        case lzl: {
          const auto& c = huff.codewords[cmd.len + lz_huff_offset];
          ret.write<bnh>({size_t(c.bitlen), c.val});
          size_t d = adr - cmd.lz_ofs - 1;
          if (cmd.type == lzs) {
            ret.write<b1, bnh>(false, {7, d});
          } else {
            ret.write<b1, bnh>(true, {10, d});
          }
        } break;
        default: assert(0);
        }
        adr += cmd.len;
      }
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
