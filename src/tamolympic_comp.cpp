#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> tamolympic_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum method { uncomp, uncompl, rle, rle0, lz2, lz };
  using tag = tag_ol<method>;

  const auto len_tab = std::to_array<vrange>({
    vrange(0x0003, 0x0003,  1, 0b1),
    vrange(0x0004, 0x0004,  2, 0b01),
    vrange(0x0005, 0x0005,  3, 0b001),
    vrange(0x0006, 0x0006,  4, 0b0001),
    vrange(0x0007, 0x0008,  6, 0b00001'0),
    vrange(0x0009, 0x0010,  9, 0b000000'000),
    vrange(0x0011, 0x0110, 14, 0b000001)         // 000001[]
  });

  const auto ofs_tab = std::to_array<vrange>({
    vrange(0x0001, 0x00ff,  9, 0b0),             // []0
    vrange(0x0100, 0x02ff, 11, 0b100),           // []1_0
    vrange(0x0300, 0x06ff, 13, 0b10100),         // []1_1_0
    vrange(0x0700, 0x0eff, 15, 0b1010100),       // []1_1_1_0
    vrange(0x0f00, 0x1eff, 17, 0b101010100),     // []1_1_1_1_0
    vrange(0x1f00, 0x3eff, 19, 0b10101010100),   // []1_1_1_1_1_0
    vrange(0x3f00, 0xbeff, 21, 0b1010101010100)  // []1_1_1_1_1_1__
  });

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost = dp[i].cost;

    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update(i, 0x14, 0x14 + 0xff, Linear<8, 20>(), {uncompl, 0, 0});

    rlen = encode::run_length(input, i, rlen);
    for (size_t li = 0; li < len_tab.size(); ++li) {
      const auto& l = len_tab[li];
      if (rlen < l.min) break;
      if (input[i] != 0) {
        dp.update(i, l.min, l.max, rlen, Constant<20>(), {rle, 0, li}, cost + l.bitlen);
      } else {
        dp.update(i, l.min, l.max, rlen, Constant<11>(), {rle0, 0, li}, cost + l.bitlen);
      }
    }

    const auto res_lz2 = lz_helper.find_best_closest(i, 0x8ff, 2);
    if (res_lz2.len == 2) {
      if ((i - res_lz2.ofs) < 0x100) {
        dp.update_lz(i, 2, 2, res_lz2, Constant<11>(), {lz2, 0, 0});
      } else {
        dp.update_lz(i, 2, 2, res_lz2, Constant<14>(), {lz2, 1, 0});
      }
    }

    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      2
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16_pre_h ret(2);
  ret.write<none>(none());

  const auto write_len = [&](const size_t li, const size_t len) -> void {
    const size_t min_len = len_tab[li].min;
    const size_t b = len_tab[li].bitlen;
    const size_t v = len_tab[li].val;
    assert(len >= min_len);
    if (b <= 9) ret.write<bnh>({b, v + (len - min_len)});
    else {
      ret.write<bnh>({b - 8, v}); ret.trim();
      ret.write<d8, none>(len - min_len, none());
    }
  };

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    const size_t li = cmd.type.li;
    const size_t oi = cmd.type.oi;
    const size_t d = adr - cmd.lz_ofs;

    switch (cmd.type.tag) {
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
      if (cmd.type.tag == rle0) ret.write<b1, d8, b1, b1>(false, 0, false, false);
      else ret.write<b1, d8, bnh, d8, b1>(false, 0, {2, 2}, input[adr], false);
      write_len(li, cmd.len);
    } break;

    case lz: {
      const auto min_dist = ofs_tab[oi].min - (oi == 0 ? 1 : 0);
      ret.write<b1, d8, b1>(false, d & 0xff, true);
      assert(d >= min_dist);
      size_t h = (d - min_dist) >> 8, t = 0;
      if (min_dist >= 0x3f00) t = h & 1, h >>= 1;
      for (size_t s = oi; s > 0; ) {
        --s;
        ret.write<bnh>({2, 2 + ((h >> s) & 1)});
      }
      ret.write<b1>(t);
      write_len(li, cmd.len);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.trim();
  assert(adr == input.size());
  assert(dp.total_cost() + 16 == ret.bit_length());
  write16(ret.out, 0, input.size());
  return ret.out;
}

} // namespace sfc_comp
