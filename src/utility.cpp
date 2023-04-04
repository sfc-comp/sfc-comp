#include "utility.hpp"

namespace sfc_comp {

void check_size(size_t input_size, size_t min_size, size_t max_size) {
  if (input_size < min_size || input_size > max_size) {
    throw std::invalid_argument(format("Input size (= 0x%X) should satisfy: 0x%X <= size <= 0x%X.",
                                input_size, min_size, max_size));
  }
}

namespace utility {

uint16_t crc16(std::span<const uint8_t> input) {
  uint16_t ret = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    ret ^= input[i];
    ret = (ret >> 8) ^ crc_table[ret & 0xff];
  }
  return ret;
}

uint16_t crc16(std::span<const uint8_t> input, size_t offset, size_t count) {
  return crc16(input.subspan(offset, count));
}

std::vector<uint16_t> u16_freq_table(std::span<const uint8_t> input, size_t k) {
  std::vector<size_t> counts(1 << 16);
  uint16_t v = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    v = input[i] << 8 | v >> 8;
    if (i > 0) counts[v] += 1;
  }
  std::vector<uint16_t> order(1 << 16);
  std::iota(order.begin(), order.end(), 0);
  std::partial_sort(order.begin(), order.begin() + k, order.end(),
    [&](const uint16_t a, const uint16_t b) { return counts[a] > counts[b]; });
  return std::vector<uint16_t>(order.begin(), order.begin() + k);
}

} // namespace utility

} // namespace sfc_comp
