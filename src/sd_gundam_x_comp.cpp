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
  struct node {
    constexpr auto operator <=> (const node& rhs) const = default;
    size_t depth;
    size_t priority;
    size_t code;
  };

  auto huff = encode::huffman(counts, true);
  auto& words = huff.words;
  auto& codewords = huff.code;
  if (words.empty()) return huffman();

  std::ranges::reverse(words);
  std::vector<size_t> depth_counts(codewords[words.front()].bitlen + 1);
  std::priority_queue<node> que;
  for (size_t i = 0; i < words.size(); ++i) {
    const size_t depth = codewords[words[i]].bitlen;
    depth_counts[depth] += 1;
    que.emplace(depth, 1, ~i);
  }

  std::vector<encode::node> nodes;
  size_t node_id = words.size();
  size_t priority = 0;
  size_t curr_depth = depth_counts.size() - 1;
  while (que.size() >= 2) {
    const auto left = que.top(); que.pop();
    const auto right = que.top(); que.pop();
    assert(left.depth == right.depth);
    if (left.depth < curr_depth) {
      curr_depth -= 1;
      assert(curr_depth > 0);
      if (depth_counts[curr_depth - 1] > 0) priority = 2 - priority;
      else priority = 0;
    }
    nodes.emplace_back(~left.code, ~right.code);
    que.emplace(left.depth - 1, priority, ~(node_id++));
  }

  const std::function<void(size_t, size_t)> traverse = [&](size_t id, size_t codeword) {
    if (id >= words.size()) {
      id -= words.size();
      traverse(nodes[id].left,  2 * codeword + 0);
      traverse(nodes[id].right, 2 * codeword + 1);
    } else {
      codewords[words[id]].val = codeword;
    }
  };
  traverse(node_id - 1, 0);
  return huffman(std::move(words), std::move(codewords), std::move(nodes));
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
        return lz_helper.find(i, max_ofs, lz_min_len);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  auto huff_bitlens = [&] {
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto h = encode::huffman(utility::freq_u8(input), true);
    for (const auto w : h.words) ret[w] = h.code[w].bitlen;
    return ret;
  }();

  const auto update_costs = [&](const huffman& huff, size_t penalty = 0) {
    const auto& codewords = huff.codewords;
    size_t longest_bit_size = 0;
    for (const auto w : huff.words) longest_bit_size = std::max<size_t>(longest_bit_size, codewords[w].bitlen);
    huff_bitlens.assign(huff_bitlens.size(), longest_bit_size + penalty);
    for (const auto w : huff.words) huff_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;
  size_t penalty = ilog2(2 * input.size() + 1) + 9;

  while (true) {
    solver<tag> dp(input.size());
    for (size_t i = input.size(); i-- > 0; ) {
      dp.update_matrix(i, ofs_tab, lz_lens,
        [&](size_t li) { return huff_bitlens[lz_lens[li] + huff_offset]; },
        [&](size_t oi) { return lz_memo[i][oi]; },
        [&](size_t oi, size_t) -> tag { return {lz, oi}; }
      );
      dp.update(i, 1, huff_bitlens[input[i]], {uncomp, 0});
    }

    std::array<size_t, 0x100 + lz_lens.size()> code_counts = {};
    {
      size_t adr = 0;
      for (const auto& cmd : dp.optimal_path()) {
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
      for (const auto& cmd : dp.optimal_path()) {
        const auto [tag, oi] = cmd.type;
        switch (tag) {
        case uncomp: {
          ret.write<bnh>(huff.codewords[input[adr]]);
        } break;
        case lz: {
          ret.write<bnh>(huff.codewords[cmd.len + huff_offset]);
          const auto& o = ofs_tab[oi];
          ret.write<bnh>({o.bitlen, o.val + ((adr - cmd.lz_ofs()) - o.min)});
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
