#pragma once

#include <cstddef>

#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <span>

namespace sfc_comp {

template <typename... Args>
std::string format(const char* format, Args... args) {
  size_t len = std::snprintf(nullptr, 0, format, args...);
  std::vector<char> buff(len + 1);
  std::snprintf(buff.data(), len + 1, format, args...);
  return std::string(buff.begin(), buff.begin() + len);
}

void check_size(size_t input_size, size_t min_size, size_t max_size);
void check_divisibility(size_t input_size, size_t divisor);

namespace utility {

uint16_t crc16(std::span<const uint8_t> input);
uint16_t crc16(std::span<const uint8_t> input, size_t offset, size_t count);

std::array<size_t, 256> freq_u8(std::span<const uint8_t> input);

std::vector<uint16_t> k_most_freq_u16(std::span<const uint8_t> input, size_t k);

template <typename T, size_t K>
requires std::integral<T>
std::array<T, K> k_most(std::span<const size_t> counts) {
  const size_t size = size_t(1) << (8 * sizeof(T));
  if (!(size == counts.size())) {
    throw std::runtime_error(format("The input size (= 0x%X) should be 0x%X.", counts.size(), size));
  }
  T order[size];
  std::iota(order, order + size, 0);
  std::partial_sort(order, order + K, order + size,
    [&](const T a, const T b) { return counts[a] > counts[b]; });
  std::array<T, K> ret;
  for (size_t i = 0; i < K; ++i) ret[i] = order[i];
  return ret;
}

} // namespace utility

template <typename T, size_t Extent, typename Func>
requires std::convertible_to<std::invoke_result_t<Func, size_t>, T>
constexpr std::array<T, Extent> create_array(Func&& func) {
  std::array<T, Extent> ret;
  for (size_t i = 0; i < Extent; ++i) ret[i] = func(i);
  return ret;
}

template <size_t Size>
constexpr auto inverse_map(std::span<const size_t> vals) {
  std::array<ptrdiff_t, Size> ret; ret.fill(-1);
  for (size_t i = 0; i < vals.size(); ++i) {
    size_t l = vals[i];
    if (l >= Size) throw std::logic_error("Too small array size.");
    if (ret[l] >= 0) throw std::logic_error("Found duplicate values.");
    ret[l] = i;
  }
  return ret;
}

inline constexpr uint64_t pext(uint64_t v, uint64_t mask) {
  uint64_t ret = 0;
  for (uint64_t t = 1; mask > 0; t <<= 1) {
    uint64_t lowest = mask & -mask;
    if (v & lowest) ret |= t;
    mask -= lowest;
  }
  return ret;
}

inline constexpr uint64_t pdep(uint64_t v, uint64_t mask) {
  uint64_t ret = 0;
  for (; mask > 0; v >>= 1) {
    uint64_t lowest = mask & -mask;
    if (v & 1) ret |= lowest;
    mask -= lowest;
  }
  return ret;
}

inline constexpr uint64_t masked_add(uint64_t src, uint64_t v, uint64_t mask) {
  return (src & ~mask) | pdep(pext(src, mask) + v, mask);
}

inline constexpr size_t ilog2(size_t n) {
  return std::bit_width(n) - 1;
}

inline constexpr uint64_t low_bits_mask(size_t bits) {
  // bits < 64
  return uint64_t(0x7fff'ffff'ffff'ffff) >> (63 - bits);
}

inline uint16_t read16(std::span<const uint8_t> input, size_t i) {
  return input[i] | (input[i + 1] << 8);
}

inline uint32_t read24(std::span<const uint8_t> input, size_t i) {
  return input[i] | (input[i + 1] << 8) | (input[i + 2] << 16);
}

inline uint32_t read32(std::span<const uint8_t> input, size_t i) {
  return read16(input, i) | read16(input, i + 2) << 16;
}

inline void write16(std::span<uint8_t> c, size_t i, uint32_t v) {
  c[i + 0] = v >> 0;
  c[i + 1] = v >> 8;
}

inline void write16b(std::span<uint8_t> c, size_t i, uint32_t v) {
  c[i + 0] = v >> 8;
  c[i + 1] = v >> 0;
}

inline void write24(std::span<uint8_t> c, size_t i, uint32_t v) {
  c[i + 0] = v >> 0;
  c[i + 1] = v >> 8;
  c[i + 2] = v >> 16;
}

inline void write24b(std::span<uint8_t> c, size_t i, uint32_t v) {
  c[i + 0] = v >> 16;
  c[i + 1] = v >> 8;
  c[i + 2] = v >> 0;
}

inline void write32(std::span<uint8_t> c, size_t i, uint32_t v) {
  c[i + 0] = v >> 0;
  c[i + 1] = v >> 8;
  c[i + 2] = v >> 16;
  c[i + 3] = v >> 24;
}

inline void write32b(std::span<uint8_t> c, size_t i, uint32_t v) {
  c[i + 0] = v >> 24;
  c[i + 1] = v >> 16;
  c[i + 2] = v >> 8;
  c[i + 3] = v >> 0;
}

inline constexpr uint16_t swap16(uint16_t x) {
  return (x >> 8) | (x << 8);
}

} // namespace sfc_comp
