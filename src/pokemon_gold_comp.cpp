#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> pokemon_gold_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x8000);

  enum tag {
    uncomp, rle, rle16, rle0, lz, lzh, lzv, lzs, lzhs, lzvs,
    uncompl, rlel, rle16l, rle0l, lzl, lzhl, lzvl, lzsl, lzhsl, lzvsl
  };

  lz_helper_kirby lz_helper(input);
  uncomp_helper u_helper(input.size(), 1);
  sssp_solver<tag> dp(input.size());

  size_t rlen = 0, rlen16 = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    u_helper.update(i, dp[i].cost);
    auto u1 = u_helper.find(i + 1, 1, 0x20);
    dp.update_u(i + 1, u1.len, uncomp, u1.cost + 1);
    auto u2 = u_helper.find(i + 1, 0x21, 0x400);
    dp.update_u(i + 1, u2.len, uncompl, u2.cost + 2);

    rlen = encode::run_length(input, i, rlen);
    dp.update(i, 1, 0x20, rlen, Constant<2>(), rle);
    dp.update(i, 0x21, 0x400, rlen, Constant<3>(), rlel);
    if (input[i] == 0) {
      dp.update(i, 1, 0x20, rlen, Constant<1>(), rle0);
      dp.update(i, 0x21, 0x400, rlen, Constant<2>(), rle0l);
    }
    rlen16 = encode::run_length16(input, i, rlen16);
    dp.update(i, 1, 0x20, rlen16, Constant<3>(), rle16);
    dp.update(i, 0x21, 0x400, rlen16, Constant<4>(), rle16l);

    auto res_lz = lz_helper.find(i, 0x8000, 3);
    dp.update_lz(i, 1, 0x20, res_lz, Constant<3>(), lz);
    dp.update_lz(i, 0x21, 0x400, res_lz, Constant<4>(), lzl);
    auto res_lzs = lz_helper.find(i, 0x80, 2);
    dp.update_lz(i, 1, 0x20, res_lzs, Constant<2>(), lzs);
    dp.update_lz(i, 0x21, 0x400, res_lzs, Constant<3>(), lzsl);

    auto res_lzh = lz_helper.find_h(i, 0x8000, 3);
    dp.update_lz(i, 1, 0x20, res_lzh, Constant<3>(), lzh);
    dp.update_lz(i, 0x21, 0x400, res_lzh, Constant<4>(), lzhl);
    auto res_lzhs = lz_helper.find_h(i, 0x80, 2);
    dp.update_lz(i, 1, 0x20, res_lzhs, Constant<2>(), lzhs);
    dp.update_lz(i, 0x21, 0x400, res_lzhs, Constant<3>(), lzhsl);

    auto res_lzv = lz_helper.find_v(i, 0x8000, 3);
    dp.update_lz(i, 1, 0x20, res_lzv, Constant<3>(), lzv);
    dp.update_lz(i, 0x21, 0x400, res_lzv, Constant<4>(), lzvl);
    auto res_lzvs = lz_helper.find_v(i, 0x80, 2);
    dp.update_lz(i, 1, 0x20, res_lzvs, Constant<2>(), lzvs);
    dp.update_lz(i, 0x21, 0x400, res_lzvs, Constant<3>(), lzvsl);

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x00 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rle: ret.write<d8, d8>(0x20 + cmd.len - 1, input[adr]); break;
    case rle16: ret.write<d8, d8, d8>(0x40 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case rle0: ret.write<d8>(0x60 + cmd.len - 1); break;
    case lzs: ret.write<d8, d8>(0x80 + cmd.len - 1, 0x80 + d - 1); break;
    case lz: ret.write<d8, d16b>(0x80 + cmd.len - 1, cmd.lz_ofs); break;
    case lzhs: ret.write<d8, d8>(0xa0 + cmd.len - 1, 0x80 + d - 1); break;
    case lzh: ret.write<d8, d16b>(0xa0 + cmd.len - 1, cmd.lz_ofs); break;
    case lzvs: ret.write<d8, d8>(0xc0 + cmd.len - 1, 0x80 + d - 1); break;
    case lzv: ret.write<d8, d16b>(0xc0 + cmd.len - 1, cmd.lz_ofs); break;

    case uncompl: ret.write<d16b, d8n>(0xe000 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rlel: ret.write<d16b, d8>(0xe400 + cmd.len - 1, input[adr]); break;
    case rle16l: ret.write<d16b, d8, d8>(0xe800 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case rle0l: ret.write<d16b>(0xec00 + cmd.len - 1); break;
    case lzsl: ret.write<d16b, d8>(0xf000 + cmd.len - 1, 0x80 + d - 1); break;
    case lzl: ret.write<d16b, d16b>(0xf000 + cmd.len - 1, cmd.lz_ofs); break;
    case lzhsl: ret.write<d16b, d8>(0xf400 + cmd.len - 1, 0x80 + d - 1); break;
    case lzhl: ret.write<d16b, d16b>(0xf400 + cmd.len - 1, cmd.lz_ofs); break;
    case lzvsl: ret.write<d16b, d8>(0xf800 + cmd.len - 1, 0x80 + d - 1); break;
    case lzvl: ret.write<d16b, d16b>(0xf800 + cmd.len - 1, cmd.lz_ofs); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.optimal_cost() == ret.size());
  assert(adr == input.size());
  ret.write<d8>(0xff);

  return ret.out;
}

} // namespace sfc_comp
