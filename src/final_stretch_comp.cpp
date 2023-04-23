#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> final_stretch_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  auto ret = lzss(
    input, 0x12, [](std::span<const uint8_t>) {},
    0x1000, 3, 0x12,
    4, true, true,
    [&](size_t, size_t o, size_t l) {
      size_t d = (o - 0x24) & 0x0fff;
      return (d & 0x0f00) << 4 | (l - 3) << 8 | (d & 0x00ff);
    }
  );
  write16(ret, 0, input.size());
  if (input.size() > 0) {
    write16(ret, 2, ret.size() - 4);
  } else {
    write16(ret, 2, 0x0001);
  }
  return ret;
}

} // namespace sfc_comp
