#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> tamolympic_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum Tag { uncomp, uncompl, rle, rle0, lz2, lz };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (!(tag == lz || tag == rle || tag == rle0)) return true;
      return li == rhs.li;
    }
    Tag tag;
    size_t oi, li;
  };

  const std::array<vrange, 7> len_tab = {
    vrange(0x0003, 0x0003,  1, 0x0001), // 1
    vrange(0x0004, 0x0004,  2, 0x0001), // 01
    vrange(0x0005, 0x0005,  3, 0x0001), // 001
    vrange(0x0006, 0x0006,  4, 0x0001), // 0001
    vrange(0x0007, 0x0008,  6, 0x0002), // 00001_
    vrange(0x0009, 0x0010,  9, 0x0000), // 000000___
    vrange(0x0011, 0x0110, 14, 0x0001)  // 000001[]
  };

  const std::array<vrange, 7> ofs_tab = {
    vrange(0x0001, 0x00ff,  9, 0x0000), // []0
    vrange(0x0100, 0x02ff, 11, 0x0004), // []1_0
    vrange(0x0300, 0x06ff, 13, 0x0014), // []1_1_0
    vrange(0x0700, 0x0eff, 15, 0x0054), // []1_1_1_0
    vrange(0x0f00, 0x1eff, 17, 0x0154), // []1_1_1_1_0
    vrange(0x1f00, 0x3eff, 19, 0x0554), // []1_1_1_1_1_0
    vrange(0x3f00, 0xbeff, 21, 0x1554)  // []1_1_1_1_1_1__
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

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
      [&](size_t oi, size_t li) -> CompType { return {lz, oi, li}; },
      2
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16_hasty ret; ret.write<d16>(0);
  ret.write<none>(none());

  const auto write_len = [&](const size_t li, const size_t len) -> void {
    const size_t min_len = len_tab[li].min;
    const size_t b = len_tab[li].bitlen;
    const size_t v = len_tab[li].val;
    assert(len >= min_len);
    if (b <= 9) ret.write<b8hn_h>({b, v + (len - min_len)});
    else {
      ret.write<b8hn_h>({b - 8, v}); ret.trim();
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
      ret.write<b1h>(true); ret.trim();
      ret.write<d8, none>(input[adr], none());
    } break;

    case uncompl: {
      ret.write<b1h, d8, b8hn_h, d8, b1h>(false, 0, {2, 2}, cmd.len - 0x14, true);
      for (size_t i = 0; i < cmd.len; ++i) ret.write<d8>(input[adr + i]);
    } break;

    case lz2: {
      ret.write<b1h, d8, b8hn_h>(false, d & 0xff, {1, 0});
      if (oi == 0) ret.write<b1h>(false);
      else ret.write<b1h, b8hn_h>(true, {3, (d - 0x100) >> 8});
    } break;

    case rle:
    case rle0: {
      if (cmd.type.tag == rle0) ret.write<b1h, d8, b1h, b1h>(false, 0, false, false);
      else ret.write<b1h, d8, b8hn_h, d8, b1h>(false, 0, {2, 2}, input[adr], false);
      write_len(li, cmd.len);
    } break;

    case lz: {
      const auto min_dist = ofs_tab[oi].min - (oi == 0 ? 1 : 0);
      ret.write<b1h, d8, b1h>(false, d & 0xff, true);
      assert(d >= min_dist);
      size_t h = (d - min_dist) >> 8, t = 0;
      if (min_dist >= 0x3f00) t = h & 1, h >>= 1;
      for (size_t s = oi; s > 0; ) {
        --s;
        ret.write<b8hn_h>({2, 2 + ((h >> s) & 1)});
      }
      ret.write<b1h>(t);
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
