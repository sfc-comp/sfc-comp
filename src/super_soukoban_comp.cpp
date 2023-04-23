#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_soukoban_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);
  std::vector<uint8_t> best;
  for (size_t comp_ty = 0; comp_ty < 4; ++comp_ty) {
    const size_t min_len = 3;
    const size_t max_len = min_len + ((0x10 << comp_ty) - 1);
    const size_t max_ofs = 0x1000 >> comp_ty;
    const size_t mask = max_ofs - 1;
    auto compressed = lzss(
      input, max_len, [] (std::span<const uint8_t>) {},
      max_ofs, min_len, max_len,
      1, true, true,
      [&](size_t, size_t o, size_t l) {
        size_t d = (o - max_len * 2) & mask;
        return ((d >> 8) << (12 + comp_ty)) | (l - min_len) << 8 | (d & 0x00ff);
      }
    );
    compressed[0] = comp_ty;
    if (best.empty() || compressed.size() < best.size()) {
      best = std::move(compressed);
    }
  }
  if (best.empty() || best.size() >= input.size() + 1) {
    best.resize(input.size() + 1);
    std::copy(input.begin(), input.end(), best.begin() + 1);
    best[0] = 0x04;
  }
  return best;
}

} // namespace sfc_comp
