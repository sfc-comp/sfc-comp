#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> dragon_knight_4_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);
  auto ret = lzss(
    input, 0, [] (std::span<const uint8_t>) {},
    0x800, 2, 0x21,
    2, false, false,
    [&] (size_t i, size_t o, size_t l) {
      size_t d = 0x0800 - (i - o);
      return (d & 0x00ff) << 8 | (l - 2) << 3 | (d & 0x0700) >> 8;
    }
  );
  write16(ret, 0, input.size());
  if (input.size() >= 2 && ret.size() >= input.size() + 2) {
    ret.resize(input.size() + 2);
    write16(ret, 0, 0x010001 - input.size());
    std::copy(input.begin(), input.end(), ret.begin() + 2);
  }
  return ret;
}

} // namespace sfc_comp
