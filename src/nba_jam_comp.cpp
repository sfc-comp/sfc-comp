#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> nba_jam_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);

  enum method { none, uncomp, lz, lz2, lzd };
  using tag = tag_ol<method>;

  static constexpr auto ulen_tab = to_vranges({
    {0x0001,  3, 0b1'00},
    {0x0004,  8, 0b1110'0000},
    {0x0014, 13, 0b11110'00000000 + 0x10},
    {0x0104, 21, 0b11111'0000000000000000 + 0x100}
  }, 0xffff);

  static constexpr size_t lz_min_len = 2;
  static constexpr auto len_tab = to_vranges({
              // 0b1 (len == 2)
    {0x0003,  2, 0b00},
    {0x0004,  7, 0b010'0000},
    {0x0014, 11, 0b011'00000000 + 0x10},
  }, 0x0103);

  static constexpr size_t lz3_min_len = len_tab.front().min;
  static constexpr size_t lz3_max_len = len_tab.back().max;
  static constexpr auto len_map = [&] {
    std::array<size_t, lz3_max_len + 1> ret = {}; ret.fill(-1);
    for (size_t li = 0; li < len_tab.size(); ++li) {
      for (size_t l = len_tab[li].min; l <= len_tab[li].max; ++l) ret[l] = li;
    }
    return ret;
  }();

  // (dist - len + 1)
  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  9, 0b0'00000000 + 1},
    {0x0100, 17, 0b1'0000000000000000 + 0x100}
  }, 0xffff);

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  non_overlapping_lz_helper nlz_helper(input);
  lz_helper lz_helper(input, true);
  solver<tag> dp0(input.size()), dp1(input.size());
  auto c0_0 = dp0.c<0>(len_tab.back().max);
  auto c8_1 = dp1.c<8>(ulen_tab.back().max);

  std::array<std::array<size_t, (lz3_max_len - lz3_min_len)>, ofs_tab.size()> res_lz_d = {};
  if (input.size() > 0) lz_helper.reset(input.size() - 1);
  for (size_t i = input.size(); i-- > 0; ) {
    if (i > 0) lz_helper.reset(i - 1);

    const auto res_lz2 = lz_helper.find(i, (ofs_tab[0].max + lz_min_len) - 1, lz_min_len);
    dp1.update(i, 2, 2, res_lz2, c0_0, 9, {lz2, 0, 0});

    // This update should be done together with the next updates.
    dp1.update_matrix(i, ofs_tab, len_tab, c0_0, 0,
      [&](size_t oi) { return nlz_helper.find_non_overlapping(i, (ofs_tab[oi].max + lz3_min_len) - 1); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    for (size_t oi = 0; oi < ofs_tab.size(); ++oi) {
      const auto& o = ofs_tab[oi];
      auto& res_lz = res_lz_d[oi];
      sliding_window_max<size_t, (lz3_max_len - lz3_min_len)> results;
      const size_t v = (o.max - o.min) + 1;
      for (size_t l = lz3_min_len + 1; l <= lz3_max_len; ++l) {
        const size_t j = l - (lz3_min_len + 1);
        const size_t dist = (o.max + l) - 1;
        if (dist <= i) {
          res_lz[j] = encode::lz_dist_r(input, i, dist, res_lz[j]);
          results.add(dist, std::min(res_lz[j], dist));
        }
        if (j >= v) results.pop(dist - v);
        if (results.empty()) break;
        const auto res = results.get();
        if (res.val < l) continue;
        const auto li = len_map[l];
        dp1.update_c(i, l, dp0[i + l].cost + len_tab[li].bitlen + o.bitlen, {lzd, oi, li}, i - res.key);
      }
    }
    dp0.update(i, ulen_tab, c8_1, 0, [&](size_t li) -> tag { return {uncomp, 0, li}; });
    dp0.update_c(i, 0, dp1[i].cost + 1, {none, 0, 0});

    c0_0.update(i); c8_1.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(2);

  size_t adr = 0;
  for (size_t curr = 0; adr < input.size(); curr ^= 1) {
    const auto cmd = (curr == 0) ? dp0[adr] : dp1[adr];
    const auto [tag, oi, li] = cmd.type;
    const size_t d = adr - cmd.lz_ofs();
    switch (tag) {
    case none: {
      ret.write<b1>(false);
    } break;
    case uncomp: {
      const auto& u = ulen_tab[li];
      ret.write<bnh>({u.bitlen, u.val + (cmd.len - u.min)});
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz2: {
      assert(d >= cmd.len);
      const size_t v = (d - cmd.len) + 1; assert(1 <= v && v <= 0xff);
      ret.write<b1, bnh>(true, {8, v});
    } break;
    case lz: case lzd: {
      const auto& l = len_tab[li]; assert(l.min <= cmd.len && cmd.len <= l.max);
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      const auto& o = ofs_tab[oi];
      assert(d >= cmd.len);
      const size_t v = (d - cmd.len) + 1; assert(o.min <= v && v <= o.max);
      ret.write<bnh>({o.bitlen, (o.val - o.min) + v});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp0.optimal_cost() + 16 == ret.bit_length());
  std::reverse(ret.out.begin() + 2, ret.out.end());
  write16(ret.out, 0, ret.size() + 1);
  ret.write<d8, d16>(8, input.size());
  return ret.out;
}

} // namespace sfc_comp
