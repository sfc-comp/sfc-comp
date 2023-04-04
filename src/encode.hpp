#pragma once

#include <cstdint>
#include <cstddef>

#include <span>

namespace sfc_comp {

namespace encode {

struct rle_data {
  size_t len;
  uint64_t v;
};

size_t run_length(std::span<const uint8_t> in, size_t adr, size_t prev_len, uint8_t delta=0);
rle_data run_length_delta(std::span<const uint8_t> in, size_t adr, rle_data prev);
rle_data run_length16_delta(std::span<const uint8_t> in, size_t adr, rle_data prev);
size_t run_length16(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t run_length24(std::span<const uint8_t> in, size_t adr, size_t prev_len);
size_t coupled8(std::span<const uint8_t> in, size_t adr, size_t max_len);
rle_data common_lo16_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len);
rle_data common_lo16(std::span<const uint8_t> in, size_t adr, size_t max_len);
rle_data common_hi16_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len);
rle_data common_hi16(std::span<const uint8_t> in, size_t adr, size_t max_len);
rle_data common_lo24_16_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len);
rle_data common_lo32_24_hint(std::span<const uint8_t> in, size_t adr, size_t prev_len);
rle_data common_hi8(std::span<const uint8_t> in, size_t adr, size_t max_len);
rle_data common_lo8(std::span<const uint8_t> in, size_t adr, size_t max_len);

} // namespace encode

} // namespace sfc_comp
