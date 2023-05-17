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

  lz_helper_kirby lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x400);
  auto c1 = dp.c<1>(0x400);

  size_t rlen = 0, rlen16 = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 0x20, c1, 1, uncomp);
    dp.update(i, 0x21, 0x400, c1, 2, uncompl);

    rlen = encode::run_length_r(input, i, rlen);
    dp.update(i, 1, 0x20, rlen, c0, 2, rle);
    dp.update(i, 0x21, 0x400, rlen, c0, 3, rlel);
    if (input[i] == 0) {
      dp.update(i, 1, 0x20, rlen, c0, 1, rle0);
      dp.update(i, 0x21, 0x400, rlen, c0, 2, rle0l);
    }
    rlen16 = encode::run_length16_r(input, i, rlen16);
    dp.update(i, 1, 0x20, rlen16, c0, 3, rle16);
    dp.update(i, 0x21, 0x400, rlen16, c0, 4, rle16l);

    const auto res_lz = lz_helper.find(i, 0x8000, 3);
    dp.update(i, 1, 0x20, res_lz, c0, 3, lz);
    dp.update(i, 0x21, 0x400, res_lz, c0, 4, lzl);
    const auto res_lzs = lz_helper.find(i, 0x80, 2);
    dp.update(i, 1, 0x20, res_lzs, c0, 2, lzs);
    dp.update(i, 0x21, 0x400, res_lzs, c0, 3, lzsl);

    const auto res_lzh = lz_helper.find_h(i, 0x8000, 3);
    dp.update(i, 1, 0x20, res_lzh, c0, 3, lzh);
    dp.update(i, 0x21, 0x400, res_lzh, c0, 4, lzhl);
    const auto res_lzhs = lz_helper.find_h(i, 0x80, 2);
    dp.update(i, 1, 0x20, res_lzhs, c0, 2, lzhs);
    dp.update(i, 0x21, 0x400, res_lzhs, c0, 3, lzhsl);

    const auto res_lzv = lz_helper.find_v(i, 0x8000, 3);
    dp.update(i, 1, 0x20, res_lzv, c0, 3, lzv);
    dp.update(i, 0x21, 0x400, res_lzv, c0, 4, lzvl);
    const auto res_lzvs = lz_helper.find_v(i, 0x80, 2);
    dp.update(i, 1, 0x20, res_lzvs, c0, 2, lzvs);
    dp.update(i, 0x21, 0x400, res_lzvs, c0, 3, lzvsl);

    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x00 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rle: ret.write<d8, d8>(0x20 + cmd.len - 1, input[adr]); break;
    case rle16: ret.write<d8, d8, d8>(0x40 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case rle0: ret.write<d8>(0x60 + cmd.len - 1); break;
    case lzs: ret.write<d8, d8>(0x80 + cmd.len - 1, 0x80 + d - 1); break;
    case lz: ret.write<d8, d16b>(0x80 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzhs: ret.write<d8, d8>(0xa0 + cmd.len - 1, 0x80 + d - 1); break;
    case lzh: ret.write<d8, d16b>(0xa0 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzvs: ret.write<d8, d8>(0xc0 + cmd.len - 1, 0x80 + d - 1); break;
    case lzv: ret.write<d8, d16b>(0xc0 + cmd.len - 1, cmd.lz_ofs()); break;

    case uncompl: ret.write<d16b, d8n>(0xe000 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rlel: ret.write<d16b, d8>(0xe400 + cmd.len - 1, input[adr]); break;
    case rle16l: ret.write<d16b, d8, d8>(0xe800 + cmd.len - 1, input[adr], input[adr + 1]); break;
    case rle0l: ret.write<d16b>(0xec00 + cmd.len - 1); break;
    case lzsl: ret.write<d16b, d8>(0xf000 + cmd.len - 1, 0x80 + d - 1); break;
    case lzl: ret.write<d16b, d16b>(0xf000 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzhsl: ret.write<d16b, d8>(0xf400 + cmd.len - 1, 0x80 + d - 1); break;
    case lzhl: ret.write<d16b, d16b>(0xf400 + cmd.len - 1, cmd.lz_ofs()); break;
    case lzvsl: ret.write<d16b, d8>(0xf800 + cmd.len - 1, 0x80 + d - 1); break;
    case lzvl: ret.write<d16b, d16b>(0xf800 + cmd.len - 1, cmd.lz_ofs()); break;
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
