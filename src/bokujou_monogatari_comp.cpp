#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> bokujou_monogatari_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  auto ret = lzss<writer_b8_l>(
    input, 0x22, [](std::span<const uint8_t>) {},
    0x800, 3, 0x22,
    4, true,
    [&](size_t, size_t o, size_t l) {
      size_t d = (o - 0x44) & 0x7ff;
      return (d & 0x00ff) | (l - 3) << 8 | (d & 0x0700) << 5;
    }
  );
  write32(ret, 0, input.size());
  return ret;
}

} // namespace sfc_comp
