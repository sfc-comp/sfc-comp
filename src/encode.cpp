#include "encode.hpp"
#include "utility.hpp"

namespace sfc_comp {

namespace encode {

size_t run_length_r(std::span<const uint8_t> in, size_t adr, size_t prev_len, uint8_t delta) {
  if (adr + 1 >= in.size()) return in.size() - std::min(in.size(), adr);
  return ((in[adr] + delta) & 0xff) == in[adr + 1] ? (prev_len + 1) : 1;
}

size_t run_length16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 2 >= in.size()) return in.size() - std::min(in.size(), adr);
  return (in[adr] == in[adr + 2]) ? (prev_len + 1) : 2;
}

size_t run_length_delta_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 1 >= in.size()) return 0;
  if (adr + 2 == in.size()) return 2;
  uint8_t diff = in[adr] + in[adr + 2] - 2 * in[adr + 1];
  return diff == 0 ? prev_len + 1 : 2;
}

size_t run_length16_delta_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 3 >= in.size()) return 0;
  if (adr + 5 >= in.size()) return 4;
  uint16_t diff = read16(in, adr) + read16(in, adr + 4) - 2 * read16(in, adr + 2);
  return diff == 0 ? prev_len + 2 : 4;
}

size_t run_length24_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 3 >= in.size()) return in.size() - std::min(in.size(), adr);
  return (in[adr] == in[adr + 3]) ? (prev_len + 1) : 3;
}

size_t paired_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 1 >= in.size()) return 0;
  return (in[adr] == in[adr + 1]) ? prev_len + 2 : 0;
}

size_t common_lo32_24_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 3 >= in.size()) return 0;
  if (adr + 7 >= in.size()) return 4;
  return read24(in, adr) == read24(in, adr + 4) ? (prev_len + 4) : 4;
}

size_t common_lo24_16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 2 >= in.size()) return 0;
  if (adr + 5 >= in.size()) return 3;
  return read16(in, adr) == read16(in, adr + 3) ? (prev_len + 3) : 3;
}

size_t common_lo16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 1 >= in.size()) return 0;
  if (adr + 3 >= in.size()) return 2;
  return (in[adr] == in[adr + 2]) ? (prev_len + 2) : 2;
}

size_t common_hi16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 1 >= in.size()) return 0;
  if (adr + 3 >= in.size()) return 2;
  return (in[adr + 1] == in[adr + 3]) ? (prev_len + 2) : 2;
}

size_t common_hi8_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 1 >= in.size()) return in.size() - std::min(in.size(), adr);
  return ((in[adr] ^ in[adr + 1]) & 0xf0) == 0 ? (prev_len + 1) : 1;
}

size_t common_lo8_r(std::span<const uint8_t> in, size_t adr, size_t prev_len) {
  if (adr + 1 >= in.size()) return in.size() - std::min(in.size(), adr);
  return ((in[adr] ^ in[adr + 1]) & 0x0f) == 0 ? (prev_len + 1) : 1;
}

size_t lz_dist_r(std::span<const uint8_t> in, size_t adr, size_t dist, size_t prev_len, uint8_t delta) {
  if (adr < dist) return 0;
  if (adr >= in.size()) return 0;
  return ((in[adr - dist] + delta) & 0xff) == in[adr] ? prev_len + 1 : 0;
}

size_t run_length(std::span<const uint8_t> in, size_t adr, size_t prev_len, uint8_t delta) {
  if (prev_len >= 2) return prev_len - 1;
  const size_t n = in.size();
  uint8_t v = (in[adr] + delta) & 0xff;
  size_t l = 1;
  for (; adr + l < n && in[adr + l] == v; ++l, v = (v + delta) & 0xff);
  return l;
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

} // namespace encode

} // namespace sfc_comp
