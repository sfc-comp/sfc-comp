#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> picross_np_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);
  auto ret = lzss(
    input, 0, [] (std::span<const uint8_t>) {},
    0x1000, 3, 0x12,
    2, true, true,
    [&] (size_t, size_t o, size_t l) {
      return (o & 0x0fff) << 4 | (l - 3);
    }
  );
  write16(ret, 0, input.size());
  return ret;
}

} // namespace sfc_comp
