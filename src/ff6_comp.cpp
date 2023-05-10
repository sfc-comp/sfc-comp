#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> ff6_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  auto ret = lzss<writer_b8_l>(
    input, 0x22, [](std::span<const uint8_t>) {},
    0x800, 3, 0x22,
    2, true,
    [&](size_t, size_t o, size_t l) {
      size_t d = (o - 0x44) & 0x7ff;
      return d | (l - 3) << 11;
    }
  );
  write16(ret, 0, ret.size());
  return ret;
}

} // namespace sfc_comp
