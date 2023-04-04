#include "huffman.hpp"

namespace sfc_comp {

namespace encode {

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

  using Code = std::pair<ptrdiff_t, size_t>;
  std::vector<Code> codes;

  std::function<void(size_t, size_t)> traverse = [&](size_t id, ptrdiff_t depth) {
    if (id < num_codes) {
      codes.push_back({depth, id});
    } else {
      id -= num_codes;
      traverse(nodes[id].left, depth + 1);
      traverse(nodes[id].right, depth + 1);
    }
  };
  traverse(node_id - 1, 0);

  // normalize the tree

  std::sort(codes.begin(), codes.end());

  size_t codeword = 0;
  ptrdiff_t depth = 0;
  std::vector<huffman_data> codewords(counts.size(), huffman_data({-1, 0}));
  std::vector<size_t> ordered_words(codes.size());

  const size_t mask_bit = (zero_to_one) ? 0 : 1;
  size_t mask = 0;
  for (size_t i = 0; i < codes.size(); ++i) {
    const auto& c = codes[i];
    while (depth < c.first) codeword <<= 1, mask = mask * 2 + mask_bit, ++depth;
    ordered_words[i] = words[c.second];
    codewords[ordered_words[i]] = {depth, codeword ^ mask};
    codeword++;
  }

  huffman_result ret;
  ret.words = std::move(ordered_words);
  ret.codewords = std::move(codewords);
  return ret;
}

} // namespace encode

} // namespace sfc_comp
