#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> rs3_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  auto ret = lzss<writer_b8_l>(
    input, 0, [](std::span<const uint8_t>) {},
    0x0fff, 3, 0x12,
    2, true,
    [&](size_t i, size_t o, size_t l) {
      size_t d = i - o;
      return (d & 0x00ff) | (l - 3) << 8 | (d & 0x0f00) << 4;
    }
  );
  write16(ret, 0, ret.size() - 2);
  return ret;
}

} // namespace sfc_comp
