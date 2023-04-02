#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> shin_megami_tensei2_comp(std::span<const uint8_t> input) {
  enum CompType {
    uncomp, rle0, rle0l, rle, lz, common16
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 32, Linear<1, 1>(), uncomp);
    rlen = encode::run_length(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, 1, 0x1f, rlen, Constant<1>(), rle0);
      dp.update(i, 0x20, 0x11f, rlen, Constant<2>(), rle0l);
    } else {
      dp.update(i, 2, 0x21, rlen, Constant<2>(), rle);
    }
    // should be called after run length functions.
    auto res_lz = lz_helper.find_best_closest(i, 0x400, 0x21);
    dp.update_lz(i, 2, 0x21, res_lz, Constant<2>(), lz);
    if (input[i] == 0) {
      auto common16_len = encode::common_lo16(input, i, 0x40).len;
      dp.update_k<2>(i, 2, 0x40, common16_len, [](size_t i) { return 1 + (i >> 1); }, common16);
    }
    lz_helper.add_element(i);
  }
  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x80 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rle0: ret.write<d8>(0xe0 + cmd.len - 1); break;
    case rle0l: ret.write<d8, d8>(0xff, cmd.len - 0x20); break;
    case rle: ret.write<d8, d8>(0xc0 + cmd.len - 2, input[adr]); break;
    case lz: ret.write<d16b>((cmd.len - 2) << 10 | (0x400 - d)); break;
    case common16: ret.write<d8, d8nk>(
      0xa0 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr + 1], 2}
    ); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d16b>(0x7fff);
  assert(dp.total_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
