#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> marvelous_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum tag {
    uncomp, rle, rle16, inc, lzl,
    uncompm, rlem, rle16m, incm, lzlm,
    uncompl, rlel, rle16l, incl, lzll,
  };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x10000);
  auto c1 = dp.c<1>(0x10000);

  size_t rlen = 0, rlen16 = 0, rleni = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 0x20, c1, 1, uncomp);
    dp.update(i, 0x21, 0x400, c1, 2, uncompm);
    dp.update(i, 0x401, 0x10000, c1, 3, uncompl);

    rlen = encode::run_length_r(input, i, rlen);
    dp.update(i, 1, 0x20, rlen, c0, 2, rle);
    dp.update(i, 0x21, 0x400, rlen, c0, 3, rlem);
    dp.update(i, 0x401, 0x10000, rlen, c0, 4, rlel);

    rlen16 = encode::run_length16_r(input, i, rlen16);
    dp.update(i, 1, 0x20, rlen16, c0, 3, rle16);
    dp.update(i, 0x21, 0x400, rlen16, c0, 4, rle16m);
    dp.update(i, 0x401, 0x10000, rlen16, c0, 5, rle16l);

    rleni = encode::run_length_r(input, i, rleni, 1);
    dp.update(i, 1, 0x20, rleni, c0, 2, inc);
    dp.update(i, 0x21, 0x400, rleni, c0, 3, incm);
    dp.update(i, 0x401, 0x10000, rleni, c0, 4, incl);

    const auto res_lzl = lz_helper.find(i, 0x10000, 3);
    dp.update(i, 1, 0x20, res_lzl, c0, 3, lzl);
    dp.update(i, 0x21, 0x400, res_lzl, c0, 4, lzlm);
    dp.update(i, 0x401, 0x10000, res_lzl, c0, 5, lzll);

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
    case lzl: ret.write<d8, d16b>(0x80 + cmd.len - 1, cmd.lz_ofs()); break;

    case uncompm: ret.write<d16b, d8n>(0xe000 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rlem: ret.write<d16b, d8>(0xe400 + cmd.len - 1, input[adr]); break;
    case rle16m: ret.write<d16b, d8, d8>(0xe800 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case incm: ret.write<d16b, d8>(0xec00 + cmd.len - 1, input[adr]); break;
    case lzlm: ret.write<d16b, d16b>(0xf000 + cmd.len - 1, cmd.lz_ofs()); break;

    case uncompl: ret.write<d8, d16b, d8n>(0xc0, cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rlel: ret.write<d8, d16b, d8>(0xc4, cmd.len - 1, input[adr]); break;
    case rle16l: ret.write<d8, d16b, d8, d8>(0xc8, cmd.len - 1, input[adr], input[adr + 1]); break;
    case incl: ret.write<d8, d16b, d8>(0xcc, cmd.len - 1, input[adr]); break;
    case lzll: ret.write<d8, d16b, d16b>(0xd0, cmd.len - 1, cmd.lz_ofs()); break;
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
