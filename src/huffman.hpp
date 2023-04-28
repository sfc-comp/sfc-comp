#pragma once

#include <cstddef>

#include <vector>
#include <queue>

#include <array>
#include <functional>

#include <span>

#include "encode.hpp"

namespace sfc_comp {

namespace encode {

struct node {
  size_t left;
  size_t right;
};

struct huffman_result {
  std::vector<size_t> words;
  std::vector<codeword> codewords;
};

huffman_result huffman(std::span<const size_t> counts, bool zero_to_one);

huffman_result length_limited_huffman(std::span<const size_t> counts, size_t limit, bool zero_to_one);

} // namespace encode

} // namespace sfc_comp
