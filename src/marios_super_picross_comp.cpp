#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> marios_super_picross_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss(
    input, 0x12, [] (std::span<const uint8_t>) {},
    0x1000, 3, 0x12,
    2, true, true,
    [&] (size_t, size_t o, size_t l) {
      size_t d = (o - 0x24) & 0xfff;
      return (d << 4) | (l - 3);
    }
  );
  write16(ret, 0, input.size());
  return ret;
}

} // namespace sfc_comp
