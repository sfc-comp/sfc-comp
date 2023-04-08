#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_jinsei_game_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  auto ret = lzss(
    input, 0x12, [] (std::span<const uint8_t>) {},
    0x1000, 3, 0x12,
    0, true, true,
    [&] (size_t, size_t o, size_t l) {
      size_t d = (o - 0x12 * 2) & 0xfff;
      return (d & 0x00ff) | (l - 3) << 8 | (d & 0x0f00) << 4;
    }
  );
  return ret;
}

} // namespace sfc_comp
