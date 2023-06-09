#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> wizardry6_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss<writer_b8_l>(
    input, 0, [](std::span<const uint8_t>) {},
    0x800, 3, 0x22,
    2, false,
    [&](size_t i, size_t o, size_t l) {
      size_t v = std::min<size_t>(i, 0x0800) - (i - o);
      return (l - 3) << 11 | v;
    }
  );
  write16(ret, 0, ret.size() - 2);
  return ret;
}

} // namespace sfc_comp
