#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_loopz_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0xffff);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  7, 0b10'00000 + 1},
    {0x0020, 10, 0b0'000000000},
    {0x0220, 16, 0b11'00000000000000},
  }, 0x421f);
  static constexpr std::array<size_t, ofs_tab.size()> ofs_prefix = {2, 1, 2};

  static constexpr auto len_tab = to_vranges({
    {0x0002,  2, 0b0'0},
    {0x0004,  4, 0b10'00},
    {0x0008,  7, 0b110'0000},
    {0x0017, 11, 0b111'00000000},
  }, 0x0116);
  static constexpr std::array<size_t, len_tab.size()> len_prefix = {1, 2, 3, 3};

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    u_helper.update(i, dp[i].cost);
    const auto u0 = u_helper.find(i + 1, 1, 1);
    dp.update_u(i + 1, u0.len, {uncomp, 0, 0}, u0.cost + 1);
    const auto u1 = u_helper.find(i + 1, 0x0f, 0x0f + 0x001f);
    dp.update_u(i + 1, u1.len, {uncomp, 0, 1}, u1.cost + 14);
    const auto u2 = u_helper.find(i + 1, 0x0f, 0x0f + 0x3fff);
    dp.update_u(i + 1, u2.len, {uncomp, 0, 2}, u2.cost + 23);
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      1
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_l ret(0x0e);
  if (input.size() > 0) ret.write<b1>(false);

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      if (cmd.type.li == 0) {
        ret.write<b1>(true);
      } else {
        ret.write<bnh, bnl>({4, 6}, {4, 0x0f});
        if (cmd.type.li == 1) {
          ret.write<b1, bnl>(true, {5, cmd.len - 0x0f});
        } else {
          ret.write<b1, bnl>(false, {14, cmd.len - 0x0f});
        }
      }
      ret.write<b8ln>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const auto wr = [&ret](const vrange& tp, size_t prefix_bits, size_t v) {
        const size_t x = tp.val + (v - tp.min);
        ret.write<bnh>({prefix_bits, x >> (tp.bitlen - prefix_bits)});
        ret.write<bnl>({tp.bitlen - prefix_bits, x});
      };
      ret.write<b1>(false);
      wr(len_tab[cmd.type.li], len_prefix[cmd.type.li], cmd.len);
      wr(ofs_tab[cmd.type.oi], ofs_prefix[cmd.type.oi], adr - cmd.lz_ofs);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  std::reverse(ret.out.begin() + 0x0e, ret.out.end());
  ret.write<d16b>(0x0f);

  write32b(ret.out, 0, 0x43724d21); // CrM!
  write16b(ret.out, 4, 0); // unknown
  write32b(ret.out, 6, input.size());
  write32b(ret.out, 10, ret.size() - 14);
  assert(adr == input.size());
  assert(dp.total_cost() + (input.size() > 0 ? 1 : 0) + (14 + 2) * 8 == ret.bit_length());

  return ret.out;
}

} // namespace sfc_comp
