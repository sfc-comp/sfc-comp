#include <cassert>

#include <stdexcept>

#include "huffman.hpp"
#include "utility.hpp"

namespace sfc_comp {

namespace encode {

namespace {

using Leaf = std::pair<ptrdiff_t, size_t>;

huffman_result create_codewords(const size_t n, std::vector<Leaf>& leaves, bool zero_to_one) {
  std::sort(leaves.begin(), leaves.end());

  size_t codeword = 0;
  ptrdiff_t depth = 0;
  std::vector<huffman_data> codewords(n, huffman_data({-1, 0}));
  std::vector<size_t> words(leaves.size());

  const size_t mask_bit = (zero_to_one) ? 0 : 1;
  size_t mask = 0;
  for (size_t i = 0; i < leaves.size(); ++i) {
    const auto& leaf = leaves[i];
    while (depth < leaf.first) codeword <<= 1, mask = mask * 2 + mask_bit, ++depth;
    words[i] = leaf.second;
    codewords[leaf.second] = {depth, codeword ^ mask};
    codeword++;
  }

  huffman_result ret;
  ret.words = std::move(words);
  ret.codewords = std::move(codewords);

  return ret;
}

} // namespace

huffman_result length_limited_huffman(std::span<const size_t> counts, const size_t limit, bool zero_to_one) {
  using Count = std::pair<ptrdiff_t, size_t>;
  std::vector<size_t> words;

  std::vector<Count> weights;
  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] == 0) continue;
    weights.push_back({counts[i], words.size()});
    words.push_back(i);
  }
  if (words.size() == 0) return huffman_result();

  if (((words.size() - 1) >> limit) != 0) {
    throw std::runtime_error(
      format("Length limited Huffman encoding failed. The number of codes (= %zu) exceeds %zu.", words.size(), size_t(1) << limit));
  }

  std::sort(weights.begin(), weights.end());
  const size_t num_codes = words.size();

  std::vector<huffman_node> nodes;
  size_t node_id = num_codes;

  std::vector<Count> items;
  for (size_t lv = limit; lv > 0; --lv) {
    std::vector<Count> merged_items;
    std::merge(weights.begin(), weights.end(), items.begin(), items.end(), std::back_inserter(merged_items));
    items.clear();
    for (size_t i = 0; i + 1 < merged_items.size(); i += 2) {
      const auto& left = merged_items[i + 0];
      const auto& right = merged_items[i + 1];
      items.emplace_back(left.first + right.first, node_id++);
      nodes.emplace_back(left.second, right.second);
    }
  }

  std::vector<size_t> depth(num_codes, 0);
  std::function<void(size_t)> traverse = [&](size_t id) {
    if (id < num_codes) {
      depth[id] += 1;
    } else {
      id -= num_codes;
      traverse(nodes[id].left);
      traverse(nodes[id].right);
    }
  };

  assert(items.size() == words.size() - 1);
  for (size_t i = 0; i < items.size(); ++i) traverse(items[i].second);

  std::vector<Leaf> leaves(words.size());
  for (size_t i = 0; i < words.size(); ++i) leaves[i] = Leaf(depth[i], words[i]);

  return create_codewords(counts.size(), leaves, zero_to_one);
}

huffman_result huffman(std::span<const size_t> counts, bool zero_to_one) {
  using Count = std::pair<ptrdiff_t, size_t>;
  std::priority_queue<Count, std::vector<Count>, std::greater<Count> > pqueue;

  std::vector<size_t> words;
  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] == 0) continue;
    pqueue.push(Count(counts[i], words.size()));
    words.push_back(i);
  }
  if (words.size() == 0) return huffman_result();

  const size_t num_codes = words.size();
  std::vector<huffman_node> nodes;

  size_t node_id = num_codes;
  while (pqueue.size() >= 2) {
    const auto left = pqueue.top(); pqueue.pop();
    const auto right = pqueue.top(); pqueue.pop();
    nodes.push_back({left.second, right.second});
    pqueue.push(Count(left.first + right.first, node_id++));
  }

  // create a huffman tree
  std::vector<Leaf> leaves;
  std::function<void(size_t, size_t)> traverse = [&](size_t id, ptrdiff_t depth) {
    if (id < num_codes) {
      leaves.emplace_back(depth, words[id]);
    } else {
      id -= num_codes;
      traverse(nodes[id].left, depth + 1);
      traverse(nodes[id].right, depth + 1);
    }
  };
  traverse(node_id - 1, 0);

  return create_codewords(counts.size(), leaves, zero_to_one);
}

} // namespace encode

} // namespace sfc_comp
