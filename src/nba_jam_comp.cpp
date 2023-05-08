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
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp0(input.size()), dp1(input.size(), -1);

  encode::lz_data res_lz2 = {0, 0};
  std::array<std::array<size_t, (lz3_max_len - lz3_min_len)>, ofs_tab.size()> res_lz_d = {};

  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost0 = dp0[i].cost;
    u_helper.update(i, cost0);

    dp1.update(i, 0, 0, Constant<1>(), {none, 0, 0}, cost0);
    for (size_t k = 0; k < ulen_tab.size(); ++k) {
      const auto u = u_helper.find(i + 1, ulen_tab[k].min, ulen_tab[k].max);
      if (u.len == u_helper.nlen) continue;
      dp1.update_u(i + 1, u.len, {uncomp, 0, k}, u.cost + ulen_tab[k].bitlen);
    }

    const auto cost1 = dp1[i].cost;
    res_lz2 = nlz_helper.find_non_overlapping(i, (ofs_tab[0].max + lz_min_len) - 1, res_lz2);
    dp0.update_lz(i, 2, 2, res_lz2, Constant<9>(), {lz2, 0, 0}, cost1);

    // This update should be done together with the next updates.
    dp0.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return nlz_helper.find_non_overlapping(i, (ofs_tab[oi].max + lz3_min_len) - 1); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      0, cost1
    );

    // [TODO] simplify ?
    for (size_t oi = 0; oi < ofs_tab.size(); ++oi) {
      const auto& o = ofs_tab[oi];
      auto& res_lz = res_lz_d[oi];
      sliding_window_max<size_t, (lz3_max_len - lz3_min_len)> results;
      const size_t v = (o.max - o.min) + 1;
      for (size_t l = lz3_min_len + 1; l <= lz3_max_len; ++l) {
        const size_t j = l - (lz3_min_len + 1);
        const size_t dist = (o.max + l) - 1;
        if (dist <= i) {
          res_lz[j] = encode::lz_dist(input, i, dist, res_lz[j]);
          results.add(dist, std::min(res_lz[j], dist));
        }
        if (j >= v) results.pop(dist - v);
        if (results.empty()) break;
        const auto res = results.get();
        if (res.val < l) continue;
        const auto li = len_map[l];
        dp0.update_lz(i, l, l, encode::lz_data(i - res.key, res.val), [&](size_t) {
          return len_tab[li].bitlen + o.bitlen;
        }, {lzd, oi, li}, cost1);
      }
    }
  }

  const auto [commands, min_cost] = [&] {
    using command_type = decltype(dp0)::vertex_type;
    std::vector<command_type> ret;
    size_t adr = input.size();
    const auto min_cost = std::min(dp0.total_cost(), dp1.total_cost());
    size_t curr = (dp0.total_cost() == min_cost) ? 0 : 1;
    while (adr > 0) {
      command_type cmd;
      if (curr == 0) cmd = dp0[adr], curr = 1, assert(cmd.len > 0);
      else cmd = dp1[adr], curr = 0;
      adr -= cmd.len;
      ret.emplace_back(cmd);
    }
    assert(adr == 0 && curr == 0);
    std::reverse(ret.begin(), ret.end());
    return std::make_pair(std::move(ret), min_cost);
  }();

  using namespace data_type;
  writer_b8_h ret(2);

  size_t adr = 0;
  for (const auto cmd : commands) {
    const size_t li = cmd.type.li;
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type.tag) {
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
      const size_t v = (d - cmd.len) + 1;
      assert(1 <= v && v <= 0xff);
      ret.write<b1, bnh>(true, {8, v});
    } break;
    case lz: case lzd: {
      const auto& l = len_tab[li];
      assert(l.min <= cmd.len && cmd.len <= l.max);
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      const auto& o = ofs_tab[cmd.type.oi];
      assert(d >= cmd.len);
      const size_t v = (d - cmd.len) + 1;
      assert(o.min <= v && v <= o.max);
      ret.write<bnh>({o.bitlen, (o.val - o.min) + v});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(min_cost + 16 == ret.bit_length());
  std::reverse(ret.out.begin() + 2, ret.out.end());
  write16(ret.out, 0, ret.size() + 1);
  ret.write<d8, d16>(8, input.size());
  return ret.out;
}

} // namespace sfc_comp
