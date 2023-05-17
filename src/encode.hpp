#pragma once

#include <cstdint>
#include <cstddef>

#include <span>

namespace sfc_comp {

namespace encode {

struct codeword {
  ptrdiff_t bitlen;
  size_t val;
};

size_t run_length_r(std::span<const uint8_t> in, size_t adr, size_t prev_len, uint8_t delta = 0);
size_t run_length_delta_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t run_length16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t run_length16_delta_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t run_length24_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t paired_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t common_lo32_24_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t common_lo24_16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t common_lo16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t common_hi16_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t common_lo8_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t common_hi8_r(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t lz_dist_r(std::span<const uint8_t> in, size_t adr, size_t dist, size_t prev_len, uint8_t delta = 0);

size_t run_length(std::span<const uint8_t> in, size_t adr, size_t prev_len, uint8_t delta = 0);
size_t run_length16(std::span<const uint8_t> in, size_t adr, size_t prev_len);

} // namespace encode

} // namespace sfc_comp
