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

  if (input.size() == 0) return std::vector<uint8_t>(2, 0);

  enum method {uncomp, lz};
  using tag = tag_o<method>;

  static constexpr size_t word_max_count = 0x180;
  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = lz_min_len + 0xff;
  static constexpr size_t huff_offset = 0x0100 - lz_min_len;
  static constexpr auto lz_lens = create_array<size_t, (lz_max_len - lz_min_len) + 1>([&](size_t i) {
    return lz_min_len + i;
  });

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  8, 0b0'0000000},
    {0x0081, 11, 0b1'0000000000 + 0x80}
  }, 0x0400);

  const auto lz_memo = [&] {
    std::vector<std::array<encode::lz_data, 2>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, ofs_tab, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find_best(i, max_ofs);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  auto huff_bitlens = [&] {
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto h = encode::huffman(utility::freq_u8(input), true);
    for (const auto w : h.words) ret[w] = h.codewords[w].bitlen;
    return ret;
  }();

  const auto update_costs = [&](const huffman& huff, size_t penalty = 0) {
    const auto& codewords = huff.codewords;
    size_t longest_bit_size = 0;
    for (const auto w : huff.words) longest_bit_size = std::max<size_t>(longest_bit_size, codewords[w].bitlen);
    huff_bitlens.assign(huff_bitlens.size(), std::max(penalty, longest_bit_size + 9));
    for (const auto w : huff.words) huff_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;
  size_t penalty = ilog2(2 * input.size() + 1) + 9;

  while (true) {
    sssp_solver<tag> dp(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, [&](size_t) { return huff_bitlens[input[i]]; }, {uncomp, 0});
      dp.update_lz_matrix(i, ofs_tab, lz_lens,
        [&](size_t oi) { return lz_memo[i][oi]; },
        [&](size_t li) { return huff_bitlens[lz_lens[li] + huff_offset]; },
        [&](size_t oi, size_t) -> tag { return {lz, oi}; }
      );
    }

    const auto commands = dp.commands();
    std::array<size_t, 0x100 + lz_lens.size()> code_counts = {};
    {
      size_t adr = 0;
      for (const auto& cmd : commands) {
        switch (cmd.type.tag) {
        case uncomp: code_counts[input[adr]] += 1; break;
        case lz: code_counts[cmd.len + huff_offset] += 1; break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());
      const size_t non_zeros = std::ranges::count_if(code_counts, [&](size_t c) { return c > 0; });
      if (non_zeros == 1) code_counts[0x100] += 1;
    }

    const auto huff = huffman_encode(code_counts);
    if (huff.words.size() > word_max_count) {
      size_t excess = huff.words.size() - word_max_count;
      for (size_t i = 0; i < huff.words.size(); ++i) {
        size_t w = huff.words[i];
        if (w < 0x100) continue;
        code_counts[w] = 0;
        --excess;
        if (excess == 0) break;
      }
      update_costs(huffman_encode(code_counts), penalty);
      penalty = penalty * 1.25 + 1;

      // this should not happen.
      if (penalty >= 1e9) throw std::runtime_error("Failed to compress the given input.");
    } else {
      using namespace data_type;
      writer_b8_h ret; ret.write<bnh>({16, input.size()});

      assert(huff.words.size() >= 2);
      ret.write<bnh>({9, huff.words.size()});
      for (const auto w : huff.words) ret.write<bnh>({9, w});

      size_t prev = 0;
      const auto write_child = [&](size_t v) {
        if (v == prev) ret.write<b1>(false);
        else ret.write<b1, bnh>(true, {10, v});
        prev = v + 1;
      };
      for (const auto& h : huff.table) {
        write_child(h.left);
        write_child(h.right);
      }

      size_t adr = 0;
      for (const auto& cmd : commands) {
        switch (cmd.type.tag) {
        case uncomp: {
          ret.write<bnh>(huff.codewords[input[adr]]);
        } break;
        case lz: {
          ret.write<bnh>(huff.codewords[cmd.len + huff_offset]);
          const auto& o = ofs_tab[cmd.type.oi];
          ret.write<bnh>({o.bitlen, o.val + ((adr - cmd.lz_ofs) - o.min)});
        } break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());

      if (best.empty() || ret.size() < best.size()) {
        best = std::move(ret.out);
        update_costs(huff);
      } else {
        break;
      }
    }
  }

  return best;
}

} // namespace sfc_comp
