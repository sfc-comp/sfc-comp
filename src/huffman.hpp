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

struct huffman_t {
  std::vector<size_t> words;
  std::vector<codeword> code;
};

huffman_t huffman(std::span<const size_t> counts, bool zero_to_one);

huffman_t length_limited_huffman(std::span<const size_t> counts, size_t limit, bool zero_to_one);

} // namespace encode

} // namespace sfc_comp
