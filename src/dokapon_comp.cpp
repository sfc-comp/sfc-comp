#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> dokapon_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x8000);

  std::vector<uint8_t> best;
  for (size_t comp_ty = 0; comp_ty < 8; ++comp_ty) {
    const size_t min_len = 2;
    const size_t max_len = min_len + ((2 << comp_ty) - 1);
    const size_t max_ofs = 0x8000 >> comp_ty;
    auto compressed = lzss(
      input, 0x0101, [] (std::span<const uint8_t>) {},
      max_ofs, min_len, max_len,
      3, true, true,
      [&] (size_t i, size_t o, size_t l) {
        size_t d = i - o;
        return ((d - 1) & 0x00ff) << 8 | (l - 2) << (7 - comp_ty) | (d - 1) >> 8;
      }
    );
    write16(compressed, 0, input.size());
    compressed[2] = comp_ty + 1;
    if (best.size() == 0 || compressed.size() < best.size()) {
      best = std::move(compressed);
    }
  }
  return best;
}

} // namespace sfc_comp
