#pragma once

#include <cstddef>

#include <vector>
#include <queue>

#include <array>
#include <functional>

#include <span>

namespace sfc_comp {

namespace encode {

struct huffman_node {
  size_t left;
  size_t right;
};

struct huffman_data {
  ptrdiff_t bit_count;
  size_t val;
};

struct huffman_result {
  std::vector<size_t> words;
  std::vector<huffman_data> codewords;
};

huffman_result huffman(std::span<const size_t> counts, bool zero_to_one);

} // namespace encode

} // namespace sfc_comp
