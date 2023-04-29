#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> tenchi_wo_kurau_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss<writer_b8_l>(
    input, 0, [] (std::span<const uint8_t>) {},
    0x0fff, 3, 0x12,
    2, true,
    [&] (size_t i, size_t o, size_t l) {
      size_t d = (i - o);
      return (l - 3) << 12 | d;
    }
  );
  write16(ret, 0, input.size());
  return ret;
}

} // namespace sfc_comp
