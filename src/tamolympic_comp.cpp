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
      return len_no == rhs.len_no;
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  struct tup {
    size_t min_len;
    size_t val;
    size_t bits;
  };

  static constexpr std::array<tup, 8>  len_tab = {
    tup({0x0003, 0x0001,  1}), // 1
    tup({0x0004, 0x0001,  2}), // 01
    tup({0x0005, 0x0001,  3}), // 001
    tup({0x0006, 0x0001,  4}), // 0001
    tup({0x0007, 0x0002,  6}), // 00001_
    tup({0x0009, 0x0000,  9}), // 000000___
    tup({0x0011, 0x0001, 14}), // 000001[]
    tup({0x0111, 0x0000,  0})
  };

  static constexpr std::array<tup, 9> ofs_tab = {
    tup({0x0001, 0x0000,  9}), // []0
    tup({0x0100, 0x0004, 11}), // []1_0
    tup({0x0300, 0x0014, 13}), // []1_1_0
    tup({0x0700, 0x0054, 15}), // []1_1_1_0
    tup({0x0f00, 0x0154, 17}), // []1_1_1_1_0
    tup({0x1f00, 0x0554, 19}), // []1_1_1_1_1_0
    tup({0x3f00, 0x1554, 21}), // []1_1_1_1_1_1__
    tup({0xbf00, 0x0000,  0}), //
  };

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost = dp[i].cost;

    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update(i, 0x14, 0x14 + 0xff, Linear<8, 20>(), {uncompl, 0, 0});

    rlen = encode::run_length(input, i, rlen);
    for (size_t li = 0; li < len_tab.size() - 1; ++li) {
      const size_t min_len = len_tab[li].min_len;
      if (min_len > rlen) break;
      const size_t max_len = len_tab[li + 1].min_len - 1;
      if (input[i] != 0) {
        dp.update(i, min_len, max_len, rlen, Constant<20>(), {rle, 0, li}, cost + len_tab[li].bits);
      } else {
        dp.update(i, min_len, max_len, rlen, Constant<11>(), {rle0, 0, li}, cost + len_tab[li].bits);
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

    for (ptrdiff_t oi = ofs_tab.size() - 2; oi >= 0; --oi) {
      const size_t min_ofs = ofs_tab[oi].min_len;
      if (min_ofs > i) continue;
      const size_t max_ofs = ofs_tab[oi + 1].min_len - 1;
      const size_t ofs_bitsize = ofs_tab[oi].bits;
      const auto res_lz = lz_helper.find_best(i, max_ofs);
      if (res_lz.len <= 2) break;
      if ((i - res_lz.ofs) < min_ofs) continue;
      for (size_t li = 0; li < len_tab.size() - 1; ++li) {
        const size_t min_len = len_tab[li].min_len;
        if (min_len > res_lz.len) break;
        const size_t max_len = len_tab[li + 1].min_len - 1;
        const size_t len_bitsize = len_tab[li].bits;
        dp.update_lz(i, min_len, max_len, res_lz,
                     Constant<2>(), {lz, size_t(oi), li}, cost + ofs_bitsize + len_bitsize);
      }
    }

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16_hasty ret; ret.write<d16>(0);
  ret.write<none>(none());

  const auto write_len = [&](const size_t li, const size_t len) -> void {
    const size_t min_len = len_tab[li].min_len;
    const size_t b = len_tab[li].bits;
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
    const size_t li = cmd.type.len_no;
    const size_t oi = cmd.type.ofs_no;
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
      const auto min_dist = ofs_tab[oi].min_len - (oi == 0 ? 1 : 0);
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
  assert((dp.total_cost() + 7) / 8 + 2 <= ret.size() && ret.size() <= (dp.total_cost() + 7) / 8 + 4);

  write16(ret.out, 0, input.size());

  return ret.out;
}

} // namespace sfc_comp
