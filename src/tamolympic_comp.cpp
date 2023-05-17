#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> tamolympic_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum method { uncomp, uncompl, rle, rle0, lz2, lz };
  using tag = tag_ol<method>;

  const auto len_tab = to_vranges({
    {0x0003,  1, 0b1},
    {0x0004,  2, 0b01},
    {0x0005,  3, 0b001},
    {0x0006,  4, 0b0001},
    {0x0007,  6, 0b00001'0},
    {0x0009,  9, 0b000000'000},
    {0x0011, 14, 0b000001'00000000}
  }, 0x0110);

  const auto ofs_tab = to_vranges({
    {0x0001,  9, 0b0,             0b0},
    {0x0100, 11, 0b100,           0b010},
    {0x0300, 13, 0b10100,         0b01010},
    {0x0700, 15, 0b1010100,       0b0101010},
    {0x0f00, 17, 0b101010100,     0b010101010},
    {0x1f00, 19, 0b10101010100,   0b01010101010},
    {0x3f00, 21, 0b1010101010100, 0b0101010101011}
  }, 0xbeff);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(len_tab.back().max);
  auto c8 = dp.c<8>(0x14 + 0xff);

  size_t rlen = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 9, {uncomp, 0, 0});
    dp.update(i, 0x14, 0x14 + 0xff, c8, 20, {uncompl, 0, 0});

    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] != 0) {
      dp.update(i, len_tab, rlen, c0, 20, [&](size_t li) -> tag { return {rle, 0, li}; });
    } else {
      dp.update(i, len_tab, rlen, c0, 11, [&](size_t li) -> tag { return {rle0, 0, li}; });
    }

    const auto res_lz2 = lz_helper.find_closest(i, 0x8ff, 2, 2);
    if (res_lz2.len == 2) {
      if ((i - res_lz2.ofs) < 0x100) {
        dp.update(i, 2, 11, {lz2, 0, 0}, res_lz2.ofs);
      } else {
        dp.update(i, 2, 14, {lz2, 1, 0}, res_lz2.ofs);
      }
    }

    dp.update_matrix(i, ofs_tab, len_tab, c0, 2,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );

    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b16_pre_h ret(2);
  ret.write<none>(none());

  const auto write_len = [&](const size_t li, const size_t len) -> void {
    const auto& l = len_tab[li];
    assert(len >= l.min);
    const size_t v = len_tab[li].val + (len - l.min);
    if (l.bitlen <= 9) ret.write<bnh>({l.bitlen, v});
    else {
      ret.write<bnh>({l.bitlen - 8, v >> 8}); ret.trim();
      ret.write<d8, none>(v & 0xff, none());
    }
  };

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    const size_t d = adr - cmd.lz_ofs();

    switch (tag) {
    case uncomp: {
      ret.write<b1>(true); ret.trim();
      ret.write<d8, none>(input[adr], none());
    } break;

    case uncompl: {
      ret.write<b1, d8, bnh, d8, b1>(false, 0, {2, 2}, cmd.len - 0x14, true);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;

    case lz2: {
      ret.write<b1, d8, bnh>(false, d & 0xff, {1, 0});
      if (oi == 0) ret.write<b1>(false);
      else ret.write<b1, bnh>(true, {3, (d - 0x100) >> 8});
    } break;

    case rle:
    case rle0: {
      if (tag == rle0) ret.write<b1, d8, b1, b1>(false, 0, false, false);
      else ret.write<b1, d8, bnh, d8, b1>(false, 0, {2, 2}, input[adr], false);
      write_len(li, cmd.len);
    } break;

    case lz: {
      const auto& o = ofs_tab[oi];
      const size_t od = (d - o.min + (oi == 0 ? 1 : 0));
      ret.write<b1, d8, b1>(false, od & 0xff, true);
      ret.write<bnh>({o.bitlen - 8, masked_add(o.val, od >> 8, o.mask)});
      write_len(li, cmd.len);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.trim();
  assert(adr == input.size());
  assert(dp.optimal_cost() + 16 == ret.bit_length());
  write16(ret.out, 0, input.size());
  return ret.out;
}

} // namespace sfc_comp
