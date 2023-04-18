#include "image.hpp"

namespace sfc_comp {

namespace image {

namespace {

inline size_t truncate(size_t s, size_t k) {
  return (s / k) * k;
}

} // namespace

std::vector<uint8_t> from_indexed256_8_1(std::span<const uint8_t> in) {
  std::vector<uint8_t> ret(truncate(in.size() / 2, 32));
  for (size_t i = 0; i < ret.size(); i += 0x20) {
    from_indexed256_8_1(in.subspan(2 * i, 0x40), std::span(ret.begin() + i, 0x20));
  }
  return ret;
}

std::vector<uint8_t> to_indexed256_8_1(std::span<const uint8_t> in) {
  std::vector<uint8_t> ret(truncate(in.size(), 32) * 2);
  for (size_t i = 0; i < ret.size() / 2; i += 0x20) {
    to_indexed256_8_1(in.subspan(i, 0x20), std::span(ret.begin() + 2 * i, 0x40));
  }
  return ret;
}

std::vector<uint8_t> to_indexed16_h_8_1(std::span<const uint8_t> in) {
  std::vector<uint8_t> ret(truncate(in.size(), 32));
  for (size_t i = 0; i < ret.size(); i += 0x20) {
    to_indexed16_h_8_1(in.subspan(i, 0x20), std::span(ret.begin() + i, 0x20));
  }
  return ret;
}

std::vector<uint8_t> to_indexed16_h_2_4(std::span<const uint8_t> in) {
  std::vector<uint8_t> ret(truncate(in.size(), 32));
  for (size_t i = 0; i < ret.size(); i += 0x20) {
    to_indexed16_h_2_4(in.subspan(i, 0x20), std::span(ret.begin() + i, 0x20));
  }
  return ret;
}

} // namespace image

} // namespace sfc_comp
