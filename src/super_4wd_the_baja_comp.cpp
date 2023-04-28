#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_4wd_the_baja_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);
  enum tag { uncomp, uncompl, lz2, lz3, lz4, lz };

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    u_helper.update(i, dp[i].cost);
    const auto u1 = u_helper.find(i + 1, 1, 8);
    dp.update_u(i + 1, u1.len, uncomp, u1.cost + 5);
    const auto u2 = u_helper.find(i + 1, 9, 9 + 0xff);
    dp.update_u(i + 1, u2.len, uncompl, u2.cost + 11);

    const auto res_lz = lz_helper.find_best(i, 0xff);
    dp.update_lz(i, 2, 2, res_lz, Constant<10>(), lz2);
    const auto res_lz3 = lz_helper.find_best(i, 0x1ff);
    dp.update_lz(i, 3, 3, res_lz3, Constant<12>(), lz3);
    const auto res_lz4 = lz_helper.find_best(i, 0x3ff);
    dp.update_lz(i, 4, 4, res_lz4, Constant<13>(), lz4);
    dp.update_lz(i, 5, 5 + 0xff, res_lz, Constant<19>(), lz);

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16 ret; ret.write<d16, b1l>(0, false);

  const auto write8 = [&ret](size_t v) {
    if (ret.bit == 0 || ret.bit >= 8) ret.write<b8ln_l>({8, v});
    else ret.write<b8ln_h>({8, v});
  };

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp:
    case uncompl: {
      if (cmd.type == uncomp) ret.write<b8ln_h, b8ln_h>({2, 0}, {3, cmd.len - 1});
      else { ret.write<b8ln_h>({3, 7}); write8(cmd.len - 9); }
      for (size_t i = 0; i < cmd.len; ++i) write8(input[adr + i]);
    } break;
    case lz2: { ret.write<b8ln_h>({2, 1}); write8(d); } break;
    case lz3: { ret.write<b8ln_h, b8ln_h>({3, 4}, {9, d}); } break;
    case lz4: { ret.write<b8ln_h, b8ln_h>({3, 5}, {10, d}); } break;
    case lz:  { ret.write<b8ln_h>({3, 6}); write8(cmd.len - 5); write8(d); } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.total_cost() + 1 + 16 == ret.bit_length());

  if (ret.size() < 4) throw std::logic_error("Input is empty.");

  write16(ret.out, 2, 0x8000 | read16(ret.out, 2) >> 1);
  assert(ret.size() % 2 == 0);
  std::reverse(ret.out.begin() + 2, ret.out.end());
  for (size_t i = 2; i < ret.size(); i += 2) std::swap(ret[i], ret[i + 1]);
  ret.write<d16>(0);
  write16(ret.out, 0, ret.size() - 4);
  write16(ret.out, ret.size() - 2, input.size());

  return ret.out;
}

} // namespace sfc_comp
