#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> fe3_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum tag {
    uncomp, rle, rle16, inc, lzl, lzs, lzcl, lzcs,
    uncompl, rlel, rle16l, incl, lzll, lzsl, lzcll, lzcsl
  };

  lz_helper_c lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x400);
  auto c1 = dp.c<1>(0x400);

  size_t rlen = 0, rlen16 = 0, rleni = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 0x20, c1, 1, uncomp);
    dp.update(i, 0x21, 0x400, c1, 2, uncompl);

    rlen = encode::run_length_r(input, i, rlen);
    dp.update(i, 1, 0x20, rlen, c0, 2, rle);
    dp.update(i, 0x21, 0x400, rlen, c0, 3, rlel);
    rlen16 = encode::run_length16_r(input, i, rlen16);
    dp.update(i, 1, 0x20, rlen16, c0, 3, rle16);
    dp.update(i, 0x21, 0x400, rlen16, c0, 4, rle16l);
    rleni = encode::run_length_r(input, i, rleni, 1);
    dp.update(i, 1, 0x20, rleni, c0, 2, inc);
    dp.update(i, 0x21, 0x400, rleni, c0, 3, incl);

    auto res_lzl = lz_helper.find(i, 0x10000, 3);
    dp.update(i, 1, 0x20, res_lzl, c0, 3, lzl);
    dp.update(i, 0x21, 0x400, res_lzl, c0, 4, lzll);
    auto res_lzcl = lz_helper.find_c(i, 0x10000, 3);
    dp.update(i, 1, 0x20, res_lzcl, c0, 3, lzcl);
    dp.update(i, 0x21, 0x400, res_lzcl, c0, 4, lzcll);
    auto res_lz = lz_helper.find(i, 0xff, 2);
    dp.update(i, 1, 0x20, res_lz, c0, 2, lzs);
    dp.update(i, 0x21, 0x400, res_lz, c0, 3, lzsl);
    auto res_lzc = lz_helper.find_c(i, 0xff, 3);
    dp.update(i, 1, 0x300, res_lzc, c0, 3, lzcsl);

    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x00 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rle: ret.write<d8, d8>(0x20 + cmd.len - 1, input[adr]); break;
    case rle16: ret.write<d8, d8, d8>(0x40 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case inc: ret.write<d8, d8>(0x60 + cmd.len - 1, input[adr]); break;
    case lzl: ret.write<d8, d16>(0x80 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzcl: ret.write<d8, d16>(0xa0 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzs: ret.write<d8, d8>(0xc0 + cmd.len - 1, adr - cmd.lz_ofs()); break;

    case uncompl: ret.write<d16b, d8n>(0xe000 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rlel: ret.write<d16b, d8>(0xe400 + cmd.len - 1, input[adr]); break;
    case rle16l: ret.write<d16b, d8, d8>(0xe800 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case incl: ret.write<d16b, d8>(0xec00 + cmd.len - 1, input[adr]); break;
    case lzll: ret.write<d16b, d16>(0xf000 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzcll: ret.write<d16b, d16>(0xf400 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzsl: ret.write<d16b, d8>(0xf800 + cmd.len - 1, adr - cmd.lz_ofs()); break;
    case lzcsl: ret.write<d16b, d8>(0xfc00 + cmd.len - 1, adr - cmd.lz_ofs()); break;
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
