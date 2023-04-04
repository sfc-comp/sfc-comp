#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> hal_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum CompType {
    uncomp, rle, rle16, inc, lz, lzh, lzv,
    uncompl, rlel, rle16l, incl, lzl, lzhl, lzvl
  };

  lz_helper_kirby lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  size_t rlen = 0, rlen16 = 0, rleni = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 0x20, Linear<1, 1>(), uncomp);
    dp.update(i, 0x21, 0x400, Linear<1, 2>(), uncompl);
    rlen = encode::run_length(input, i, rlen);
    dp.update(i, 1, 0x20, rlen, Constant<2>(), rle);
    dp.update(i, 0x21, 0x400, rlen, Constant<3>(), rlel);
    rlen16 = encode::run_length16(input, i, rlen16);
    dp.update_k<2>(i, 2, 0x40, rlen16, Constant<3>(), rle16);
    dp.update_k<2>(i, 0x42, 0x800, rlen16, Constant<4>(), rle16l);
    rleni = encode::run_length(input, i, rleni, 1);
    dp.update(i, 1, 0x20, rleni, Constant<2>(), inc);
    dp.update(i, 0x21, 0x400, rleni, Constant<3>(), incl);

    auto res_lz = lz_helper.find_best(i, 0x10000);
    dp.update_lz(i, 1, 0x20, res_lz, Constant<3>(), lz);
    dp.update_lz(i, 0x21, 0x400, res_lz, Constant<4>(), lzl);
    auto res_lzh = lz_helper.find_best_h(i, 0x10000);
    dp.update_lz(i, 1, 0x20, res_lzh, Constant<3>(), lzh);
    dp.update_lz(i, 0x21, 0x400, res_lzh, Constant<4>(), lzhl);
    auto res_lzv = lz_helper.find_best_v(i, 0x10000);
    dp.update_lz(i, 1, 0x20, res_lzv, Constant<3>(), lzv);
    dp.update_lz(i, 0x21, 0x400, res_lzv, Constant<4>(), lzvl);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x00 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rle: ret.write<d8, d8>(0x20 + cmd.len - 1, input[adr]); break;
    case rle16: ret.write<d8, d8, d8>(0x40 + (cmd.len >> 1) - 1, input[adr], input[adr + 1]); break;
    case inc: ret.write<d8, d8>(0x60 + cmd.len - 1, input[adr]); break;
    case lz: ret.write<d8, d16>(0x80 + cmd.len - 1, cmd.lz_ofs); break;
    case lzh: ret.write<d8, d16>(0xa0 + cmd.len - 1, cmd.lz_ofs); break;
    case lzv: ret.write<d8, d16>(0xc0 + cmd.len - 1, cmd.lz_ofs); break;

    case uncompl: ret.write<d16b, d8n>(0xe000 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rlel: ret.write<d16b, d8>(0xe400 + cmd.len - 1, input[adr]); break;
    case rle16l: ret.write<d16b, d8, d8>(0xe800 + (cmd.len >> 1) - 1, input[adr], input[adr + 1]); break;
    case incl: ret.write<d16b, d8>(0xec00 + cmd.len - 1, input[adr]); break;
    case lzl: ret.write<d16b, d16>(0xf000 + cmd.len - 1, cmd.lz_ofs); break;
    case lzhl: ret.write<d16b, d16>(0xf400 + cmd.len - 1, cmd.lz_ofs); break;
    case lzvl: ret.write<d16b, d16>(0xf800 + cmd.len - 1, cmd.lz_ofs); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.total_cost() == ret.size());
  assert(adr == input.size());
  ret.write<d8>(0xff);

  return ret.out;
}

} // namespace sfc_comp
