#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> gokinjo_boukentai_comp_1(std::span<const uint8_t> input, const size_t header_size) {
  check_size(input.size(), 1, 0xffff);

  const auto huff = encode::huffman(utility::freq_u8(input), true);
  using namespace data_type;
  writer_b8_h ret(header_size);

  for (const auto w : huff.words) {
    const auto c = huff.codewords[w];
    const size_t r_zero = std::min<size_t>(std::countr_zero(c.val), c.bitlen);
    ret.write<bnh, b1, bnh>({r_zero, low_bits_mask(r_zero)}, false, {8, w});
  }
  for (const auto v : input) ret.write<bnh>(huff.codewords[v]);
  return ret.out;
}

std::vector<uint8_t> gokinjo_boukentai_comp_2(std::span<const uint8_t> input, const size_t header_size) {
  check_size(input.size(), 0, 0x800000);

  using namespace data_type;
  writer ret;
  for (size_t i = 0; i < header_size; ++i) ret.write<d8>(0);

  size_t last_len = 0;
  for (size_t i = 0; i < input.size(); ) {
    const auto v = input[i];
    size_t len = encode::run_length(input, i, 0);
    i += len;
    while (len >= 2) {
      const size_t l = std::min<size_t>(len, 256);
      ret.write<d8, d8, d8>(v, v, l - 1);
      len -= l;
    }
    last_len = len;
    if (len == 1) ret.write<d8>(v);
  }
  const auto terminator = (last_len == 1 && input.back() == 0) ? 1 : 0;
  ret.write<d8, d8, d8>(terminator, terminator, 0);
  return ret.out;
}

std::vector<uint8_t> gokinjo_boukentai_comp_3(std::span<const uint8_t> input, const size_t header_size) {
  check_size(input.size(), 0, 0xffff);

  using namespace data_type;
  writer ret;
  for (size_t i = 0; i < header_size; ++i) ret.write<d8>(0);

  size_t last_len = 0;
  for (size_t i = 0; i + 1 < input.size(); ) {
    const auto v = read16(input, i);
    size_t len = encode::run_length16(input, i, 0) / 2 * 2;
    i += len;
    while (len >= 4) {
      const size_t l = std::min<size_t>(len, 512);
      ret.write<d16, d16, d8>(v, v, l / 2 - 1);
      len -= l;
    }
    last_len = len;
    if (len == 2) ret.write<d16>(v);
  }
  if (input.size() & 1) {
    uint16_t v = input.back();
    if (last_len == 2) {
      // e.g.
      // - 0x01 0x00 0x01
      // - 0x00 0x01 0x00
      const uint16_t last_val = read16(input, input.size() - 3);
      while (v == last_val || v == 0) v += 0x0100;
      ret.write<d16, d16>(v, 0);
    } else {
      // e.g.
      // - 0x00
      // - 0x01
      // - 0x00 0x01 0x00 0x01 0x00
      if (v == 0) v += 0x100;
      ret.write<d16, d16>(v, 0);
    }
  } else {
    // The original compressor doesn't seem to consider this case,
    // or it assumes that the first two bytes of the file header (= 0x4159) would work as a terminator.
    if (last_len == 2) {
      const uint16_t last_val = read16(input, input.size() - 2);
      uint16_t v = 0;
      while (v == last_val) v += 1;
      ret.write<d16>(v);
    }
  }
  return ret.out;
}

constexpr auto file_header = std::to_array<uint8_t>({0x59, 0x41, 0x4D, 0x43, 0x52, 0x4D}); // YAMCRM

} // namespace

std::vector<uint8_t> gokinjo_boukentai_comp(std::span<const uint8_t> input) {
  static constexpr size_t header_size = 16;

  size_t comp_type = 3;
  auto best = gokinjo_boukentai_comp_3(input, header_size);

  if (auto res = gokinjo_boukentai_comp_2(input, header_size); res.size() < best.size()) {
    best = std::move(res); comp_type = 2;
  }

  if (input.size() > 0) {
    if (auto res = gokinjo_boukentai_comp_1(input, header_size); res.size() < best.size()) {
      best = std::move(res); comp_type = 1;
    }
  }

  std::copy(file_header.begin(), file_header.end(), best.begin());
  best[6] = comp_type;
  write16(best, 7, input.size());

  return best;
}

std::vector<uint8_t> bounty_sword_comp(std::span<const uint8_t> input) {
  static constexpr size_t header_size = 20;
  auto best = gokinjo_boukentai_comp_1(input, header_size);
  std::copy(file_header.begin(), file_header.end(), best.begin());
  write16(best, 16, input.size());
  return best;
}

} // namespace sfc_comp
