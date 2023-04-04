#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> dq6_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  auto ret = lzss(
    input, 0x42, [] (std::span<const uint8_t>) {},
    0x0400, 3, 0x42,
    0, true, true,
    [&] (size_t, size_t o, size_t l) {
      size_t d = (o - 0x42 * 2) & 0x03ff;
      return (d & 0x00ff) | (l - 3) << 8 | ((d & 0x0300) << 6);
    }
  );
  return ret;
}

} // namespace sfc_comp
