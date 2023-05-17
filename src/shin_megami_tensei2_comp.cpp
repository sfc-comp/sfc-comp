#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> shin_megami_tensei2_comp(std::span<const uint8_t> input) {
  enum tag { uncomp, rle0, rle0l, rle, lz, common16 };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x11f);
  auto c1 = dp.c<1>(32);
  auto c1_2 = dp.c<1, 2>(0x40);

  size_t rlen = 0; std::array<size_t, 2> lo16_lens = {};
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 32, c1, 1, uncomp);
    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, 1, 0x1f, rlen, c0, 1, rle0);
      dp.update(i, 0x20, 0x11f, rlen, c0, 2, rle0l);
    } else {
      dp.update(i, 2, 0x21, rlen, c0, 2, rle);
    }
    const auto res_lz = lz_helper.find(i, 0x400, 2);
    if (res_lz.len > 0 && (i - res_lz.ofs) >= 2) {
      dp.update(i, 2, 0x21, res_lz, c0, 2, lz);
    }
    lo16_lens[i & 1] = encode::common_lo16_r(input, i, lo16_lens[i & 1]);
    if (input[i] == 0) dp.update(i, 2, 0x40, lo16_lens[i & 1], c1_2, 1, common16);
    c0.update(i); c1.update(i); c1_2.update(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case lz: {
      assert(!(d == 1 && cmd.len == 0x21)); // 0x7fff
      ret.write<d16b>((cmd.len - 2) << 10 | (0x400 - d));
    } break;
    case uncomp: ret.write<d8, d8n>(0x80 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case common16: ret.write<d8, d8nk>(0xa0 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr + 1], 2}); break;
    case rle: ret.write<d8, d8>(0xc0 + cmd.len - 2, input[adr]); break;
    case rle0: ret.write<d8>(0xe0 + cmd.len - 1); break;
    case rle0l: ret.write<d8, d8>(0xff, cmd.len - 0x20); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d16b>(0x7fff);
  assert(dp.optimal_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
