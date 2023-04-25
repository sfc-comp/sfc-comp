#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_loopz_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0xffff);

  enum Tag {
    uncomp, lz
  };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return len_no == rhs.len_no;
      return len_no == rhs.len_no;
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  struct tup {
    size_t min_len;
    size_t val;
    size_t prefix_bits;
    size_t bits;
  };

  static constexpr std::array<tup, 5> len_tab = {
    tup({0x0002, 0x0000,  1, 2}),  // 0_
    tup({0x0004, 0x0008,  2, 4}),  // 10__
    tup({0x0008, 0x0060,  3, 7}),  // 110____
    tup({0x0017, 0x0700,  3, 11}), // 111________
    tup({0x0117, 0x0000,  0, 0})
  };

  static constexpr std::array<tup, 4> ofs_tab = {
    tup({0x0001, 0x0041, 2, 7}),  // 10_____
    tup({0x0020, 0x0000, 1, 10}), // 0_________
    tup({0x0220, 0xc000, 2, 16}), // 11______________
    tup({0x4220, 0x0000, 0, 0})
  };

  std::vector<uint8_t> input(in.begin(), in.end());
  std::reverse(input.begin(), input.end());

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    u_helper.update(i, dp[i].cost);
    const auto u0 = u_helper.find(i + 1, 1, 1);
    dp.update_u(i + 1, u0.len, {uncomp, 0, 0}, u0.cost + 1);
    const auto u1 = u_helper.find(i + 1, 0x0f, 0x0f + 0x001f);
    dp.update_u(i + 1, u1.len, {uncomp, 0, 1}, u1.cost + 14);
    const auto u2 = u_helper.find(i + 1, 0x0f, 0x0f + 0x3fff);
    dp.update_u(i + 1, u2.len, {uncomp, 0, 2}, u2.cost + 23);

    const auto cost = dp[i].cost;
    for (ptrdiff_t oi = ofs_tab.size() - 2; oi >= 0; --oi) {
      const size_t min_ofs = ofs_tab[oi].min_len;
      if (min_ofs > i) continue;
      const size_t max_ofs = ofs_tab[oi + 1].min_len - 1;
      const size_t ofs_bitsize = ofs_tab[oi].bits;
      auto res_lz = lz_helper.find_best(i, max_ofs);
      if (res_lz.len <= 1) break;
      if ((i - res_lz.ofs) < min_ofs) continue;
      for (size_t li = 0; li < len_tab.size() - 1; ++li) {
        const size_t min_len = len_tab[li].min_len;
        if (min_len > res_lz.len) break;
        const size_t max_len = len_tab[li + 1].min_len - 1;
        const size_t len_bitsize = len_tab[li].bits;
        dp.update_lz(i, min_len, max_len, res_lz,
                     Constant<0>(), {lz, size_t(oi), li}, cost + 1 + ofs_bitsize + len_bitsize);
      }
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  for (size_t i = 0; i < 0x0e; ++i) ret.write<d8>(0);
  if (input.size() > 0) ret.write<b1l>(false);

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      if (cmd.type.len_no == 0) {
        ret.write<b1l>(true);
      } else {
        ret.write<b8ln_h, b8ln_l>({4, 6}, {4, 0x0f});
        if (cmd.type.len_no == 1) {
          ret.write<b1l, b8ln_l>(true, {5, cmd.len - 0x0f});
        } else {
          ret.write<b1l, b8ln_l>(false, {14, cmd.len - 0x0f});
        }
      }
      for (size_t i = 0; i < cmd.len; ++i) ret.write<b8ln_l>({8, input[adr + i]});
    } break;
    case lz: {
      const auto wr = [&ret](const tup& tp, size_t v) {
        const size_t x = tp.val + (v - tp.min_len);
        ret.write<b8ln_h>({tp.prefix_bits, x >> (tp.bits - tp.prefix_bits)});
        ret.write<b8ln_l>({tp.bits - tp.prefix_bits, x});
      };
      ret.write<b1h>(false);
      wr(len_tab[cmd.type.len_no], cmd.len);
      wr(ofs_tab[cmd.type.ofs_no], adr - cmd.lz_ofs);
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
  assert((dp.total_cost() + (input.size() > 0 ? 1 : 0) + 7) / 8 + 14 + 2 == ret.size());

  return ret.out;
}

} // namespace sfc_comp
