#include "algorithm.hpp"
#include "encode.hpp"
#include "writer.hpp"
#include "utility.hpp"

namespace sfc_comp {

std::vector<uint8_t> front_mission_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum tag { uncomp, rle, lz };

  lz_helper lz_helper(input, true);
  if (input.size() > 0) lz_helper.reset(input.size() - 1);

  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x82);
  auto c0_2 = dp.c<0, 2>(0x12);

  size_t res_lz1 = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    if (i > 0) lz_helper.reset(i - 1);
    dp.update(i, 1, 9, uncomp);
    res_lz1 = encode::lz_dist_r(input, i, 1, res_lz1);
    dp.update(i, 3, 0x82, res_lz1, c0, 9, rle);
    dp.update(i, 4, 0x12, lz_helper.find(i, 0xfff, 4), c0_2, 17, lz);
    c0.update(i); c0_2.update(i);
  }

  using namespace data_type;
  writer_b8_l ret(2);
  size_t adr = 0, num_codes = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: ret.write<b1, d8>(false, input[adr]); break;
    case rle: ret.write<b1, d8>(true, (cmd.len - 3) << 1); break;
    case lz: {
      assert(d >= 2);
      ret.write<b1, d16b>(true, (d & 0xf00) << 4 | (cmd.len - 4) << 8 | 0x100 | (d & 0x0ff));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
    num_codes += 1;
  }
  write16(ret.out, 0, num_codes);
  assert(adr == input.size());
  assert(dp.optimal_cost() + 2 * 8 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
