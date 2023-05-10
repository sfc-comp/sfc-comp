#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> hanjuku_hero_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x800000);
  auto ret = lzss<writer_b8_l>(
    input, 0xfee, [](std::span<uint8_t> in) {
      size_t index = 0;
      for (size_t i = 0; i < 0x0100; ++i) {
        for (size_t j = 0; j < 0x0d; ++j) in[index++] = i;
      }
      for (size_t i = 0; i < 0x100; ++i) in[index++] = i;
      for (size_t i = 0; i < 0x100; ++i) in[index++] = ~i;
      for (size_t i = 0; i < 0x80; ++i) in[index++] = 0;
      for (size_t i = 0; i < 0x6e; ++i) in[index++] = 0x20;
    },
    0x1000, 3, 0x12,
    8, true,
    [&](size_t, size_t o, size_t l) {
      return (o & 0x00ff) | (l - 3) << 8 | (o & 0x0f00) << 4;
    }
  );
  write32(ret, 0, ret.size() - 8);
  write32(ret, 4, input.size());
  return ret;
}

} // namespace sfc_comp
