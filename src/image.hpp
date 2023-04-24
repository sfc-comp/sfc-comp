#pragma once

#include <cstdint>
#include <cstdio>

#include <vector>

#include <span>

namespace sfc_comp {

namespace snes4bpp {

std::vector<uint8_t> to_indexed256_8_1(std::span<const uint8_t> in);
std::vector<uint8_t> to_indexed16_h_8_1(std::span<const uint8_t> in);
std::vector<uint8_t> to_indexed16_h_2_8(std::span<const uint8_t> in);

std::vector<uint8_t> from_indexed256_8_1(std::span<const uint8_t> in);

inline void from_indexed256_8_1(std::span<const uint8_t> in, std::span<uint8_t> out) {
  for (size_t row = 0; row < 8; ++row) {
    uint64_t v = in[8 * row + 7] <<  0 | in[8 * row + 6] <<  8
               | in[8 * row + 5] << 16 | in[8 * row + 4] << 24;
    v &= 0x0f0f0f0f;
    v |= (in[8 * row + 3] << 4  | in[8 * row + 2] << 12
         |in[8 * row + 1] << 20 | in[8 * row + 0] << 28) & 0xf0f0f0f0;
    uint64_t t = (v ^ (v >> 7))  & 0x00aa00aa; v ^= t | (t << 7);
    t = (v ^ (v >> 14)) & 0x0000cccc; v ^= t | (t << 14);
    out[2 * row + 0x00] = v >>  0; out[2 * row + 0x01] = v >> 8;
    out[2 * row + 0x10] = v >> 16; out[2 * row + 0x11] = v >> 24;
  }
}

inline void to_indexed256_8_1(std::span<const uint8_t> in, std::span<uint8_t> out) {
  for (size_t row = 0; row < 8; ++row) {
    uint64_t v = in[2 * row + 0x11] <<  0 | in[2 * row + 0x10] <<  8
               | in[2 * row + 0x01] << 16 | in[2 * row + 0x00] << 24;
    v <<= 32;
    uint64_t t = (v ^ (v >> 9)) & 0x0055005500000000; v ^= t | (t << 9);
    t = (v ^ (v >> 18)) & 0x0000333300000000; v ^= t | (t << 18);
    v = (v & 0xf0f0f0f000000000) >> 36 | (v & 0x0f0f0f0f00000000);
    out[8 * row + 0] = v >>  0; out[8 * row + 1] = v >>  8;
    out[8 * row + 2] = v >> 16; out[8 * row + 3] = v >> 24;
    out[8 * row + 4] = v >> 32; out[8 * row + 5] = v >> 40;
    out[8 * row + 6] = v >> 48; out[8 * row + 7] = v >> 56;
  }
}

inline void to_indexed16_h_8_1(std::span<const uint8_t> in, std::span<uint8_t> out) {
  for (size_t row = 0; row < 8; ++row) {
    uint32_t v = in[2 * row + 0x11] << 24 | in[2 * row + 0x10] << 16
               | in[2 * row + 0x01] <<  8 | in[2 * row + 0x00] <<  0;
    uint32_t t = (v ^ (v >> 7)) & 0x00aa00aa; v ^= t | (t << 7);
    t = (v ^ (v >> 14)) & 0x0000cccc; v ^= t | (t << 14);
    t = (v ^ (v >>  4)) & 0x00f000f0; v ^= t | (t <<  4);
    out[4 * row + 0] = v >> 24;
    out[4 * row + 1] = v >>  8;
    out[4 * row + 2] = v >> 16;
    out[4 * row + 3] = v >>  0;
  }
}

inline void to_indexed16_h_2_8(std::span<const uint8_t> in, std::span<uint8_t> out) {
  for (size_t row = 0; row < 8; ++row) {
    uint32_t v = in[2 * row + 0x11] << 24 | in[2 * row + 0x10] << 16
               | in[2 * row + 0x01] <<  8 | in[2 * row + 0x00] <<  0;
    uint32_t t = (v ^ (v >> 7)) & 0x00aa00aa; v ^= t | (t << 7);
    t = (v ^ (v >> 14)) & 0x0000cccc; v ^= t | (t << 14);
    t = (v ^ (v >>  4)) & 0x00f000f0; v ^= t | (t <<  4);
    out[row + 0x00] = v >> 24;
    out[row + 0x08] = v >>  8;
    out[row + 0x10] = v >> 16;
    out[row + 0x18] = v >>  0;
  }
}

} // namespace snes4bpp

} // namespace sfc_comp
