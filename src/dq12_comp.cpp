#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> dq12_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss(
    input, 0x12, [] (std::span<const uint8_t>) {},
    0x1000, 3, 0x12,
    2, true, true,
    [&] (size_t, size_t o, size_t l) {
      size_t d = (o - 0x24) & 0xfff;
      return (d & 0x00ff) | (l - 3) << 8 | (d & 0x0f00) << 4;
    }
  );
  write16b(ret, 0, input.size());
  return ret;
}

std::vector<uint8_t> dq5_comp_2(std::span<const uint8_t> input) {
  auto ret = dq12_comp(input);
  std::swap(ret[0], ret[1]);
  return ret;
}

} // namespace sfc_comp
