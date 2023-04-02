#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> wizardry5_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss(
    input, 0, [] (std::span<const uint8_t>) {},
    0x400, 3, 0x42,
    0, true, false,
    [&] (size_t i, size_t o, size_t l) {
      size_t v = std::min<size_t>(i, 0x0400) - (i - o);
      return (l - 3) << 10 | v;
    }
  );
  return ret;
}

} // namespace sfc_comp
