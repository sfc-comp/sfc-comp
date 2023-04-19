#include "encode.hpp"
#include "utility.hpp"

namespace sfc_comp {

namespace encode {

size_t run_length(std::span<const uint8_t> in, size_t adr, size_t prev_len, uint8_t delta) {
  if (prev_len >= 2) return prev_len - 1;
  const size_t n = in.size();
  uint8_t v = (in[adr] + delta) & 0xff;
  size_t l = 1;
  for (; adr + l < n && in[adr + l] == v; ++l, v = (v + delta) & 0xff);
  return l;
}

rle_data run_length_delta(std::span<const uint8_t> in, size_t adr, rle_data prev) {
  if (prev.len >= 3) return {prev.len - 1, prev.v};
  const size_t n = in.size();
  if (adr + 1 >= n) return {0, 0};
  uint8_t v = in[adr + 1];
  const uint8_t delta = v - in[adr];
  size_t l = 2;
  for (v += delta; adr + l < n && in[adr + l] == v; ++l, v += delta);
  return {l, delta};
}

rle_data run_length16_delta(std::span<const uint8_t> in, size_t adr, rle_data prev) {
  if (prev.len >= 6) return {prev.len - 2, prev.v};
  const size_t n = in.size();
  if (adr + 3 >= n) return {0, 0};
  uint16_t v = read16(in, adr + 2);
  const uint16_t delta = v - read16(in, adr);
  size_t l = 4;
  for (v += delta; adr + l + 1 < n && read16(in, adr + l) == v; l += 2, v += delta);
  return {l, delta};
}

size_t run_length16(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (prev_len >= 3) return prev_len - 1;
  const size_t n = in.size();
  if (adr + 1 >= n) return 0;
  uint8_t v1 = in[adr], v2 = in[adr + 1];
  size_t l = 2;
  for (; adr + l < n && in[adr + l] == v1; ++l, std::swap(v1, v2));
  return l;
}

size_t run_length24(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (prev_len >= 4) return prev_len - 1;
  const size_t n = in.size();
  if (adr + 2 >= n) return 0;
  size_t l = 3;
  for (; adr + l < n && in[adr + l] == in[adr + l - 3]; ++l);
  return l;
}

size_t coupled8(std::span<const uint8_t> in, size_t adr, size_t max_len) {
  const size_t n = in.size();
  size_t l = 0;
  for (; l + 2 <= max_len && adr + l + 1 < n && in[adr + l] == in[adr + l + 1]; l += 2);
  return l;
}

rle_data common_lo16_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (prev_len >= 4) return {prev_len - 2, in[adr]};
  const size_t n = in.size();
  if (adr + 1 >= n) return {0, 0};
  uint8_t v = in[adr];
  size_t l = 2;
  for (; adr + l + 1 < n && in[adr + l] == v; l += 2);
  return {l, v};
}

rle_data common_lo16(std::span<const uint8_t> in, size_t adr, size_t max_len) {
  const size_t n = in.size();
  if (adr + 1 >= n) return {0, 0};
  uint8_t v = in[adr];
  size_t l = 2;
  for (; l + 2 <= max_len && adr + l + 1 < n && in[adr + l] == v; l += 2);
  return {l, v};
}

rle_data common_hi16_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (prev_len >= 4) return {prev_len - 2, in[adr + 1]};
  const size_t n = in.size();
  if (adr + 1 >= n) return {0, 0};
  uint8_t v = in[adr + 1];
  size_t l = 2;
  for (; adr + l + 1 < n && in[adr + l + 1] == v; l += 2);
  return {l, v};
}

rle_data common_hi16(std::span<const uint8_t> in, size_t adr, size_t max_len) {
  const size_t n = in.size();
  if (adr + 1 >= n) return {0, 0};
  uint8_t v = in[adr + 1];
  size_t l = 2;
  for (; l + 2 <= max_len && adr + l + 1 < n && in[adr + l + 1] == v; l += 2);
  return {l, v};
}

rle_data common_lo24_16_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (prev_len >= 6) return {prev_len - 3, static_cast<size_t>(in[adr] | in[adr + 1] << 8)};
  const size_t n = in.size();
  if (adr + 2 >= n) return {0, 0};
  uint8_t v1 = in[adr], v2 = in[adr + 1];
  size_t l = 3;
  for (; adr + l + 2 < n &&
         in[adr + l] == v1 && in[adr + l + 1] == v2; l += 3);
  return {l, static_cast<size_t>(v1 | (v2 << 8))};
}

rle_data common_lo32_24_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (prev_len >= 8) return {prev_len - 4,
                             static_cast<size_t>(in[adr] | in[adr + 1] << 8 | in[adr + 2] << 16)};
  const size_t n = in.size();
  if (adr + 3 >= n) return {0, 0};
  uint8_t v1 = in[adr], v2 = in[adr + 1], v3 = in[adr + 2];
  size_t l = 4;
  for (; adr + l + 3 < n &&
         in[adr + l] == v1 && in[adr + l + 1] == v2 && in[adr + l + 2] == v3; l += 4);
  return {l, static_cast<size_t>(v1 | v2 << 8 | v3 << 16)};
}

rle_data common_hi8(std::span<const uint8_t> in, size_t adr, size_t max_len) {
  const size_t n = in.size();
  size_t v = in[adr] & 0xf0;
  size_t l = 1;
  for (; l + 1 <= max_len && adr + l < n && (in[adr + l] & 0xf0) == v; l += 1);
  return {l, v};
}

rle_data common_lo8(std::span<const uint8_t> in, size_t adr, size_t max_len) {
  const size_t n = in.size();
  size_t v = in[adr] & 0x0f;
  size_t l = 1;
  for (; l + 1 <= max_len && adr + l < n && (in[adr + l] & 0x0f) == v; l += 1);
  return {l, v};
}

size_t lz_dist(std::span<const uint8_t> in, size_t adr, size_t dist, size_t prev_len) {
  if (prev_len >= 1) return prev_len - 1;
  if (adr < dist) return 0;
  const size_t n = in.size();
  size_t l = 0;
  for (; adr + l < n && in[adr + l - dist] == in[adr + l]; l += 1);
  return l;
}

} // namespace encode

} // namespace sfc_comp
