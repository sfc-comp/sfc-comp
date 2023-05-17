#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sailor_moon_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum tag { uncomp, lzs, lzls, lzll };
  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.c<0>(256);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    dp.update(i, 2, 5, lz_helper.find(i, 0x100, 2), c0, 12, lzs);
    const auto res_lzl = lz_helper.find(i, 0x2000, 3);
    dp.update(i, 3, 9, res_lzl, c0, 18, lzls);
    dp.update(i, 10, 256, res_lzl, c0, 26, lzll);
    c0.update(i);
  }

  using namespace data_type;
  writer_b16_pre_l ret;
  ret.write<none>(none());

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: {
      ret.write<bnh, d8>({1, 1}, input[adr]);
    } break;
    case lzs: {
      ret.write<bnh, d8>({4, cmd.len - 2}, 0x100 - d);
    } break;
    case lzls: case lzll: {
      ret.write<bnh, d8>({2, 1}, (0x2000 - d) & 0x00ff);
      if (cmd.type == lzls) {
        ret.write<d8>((cmd.len - 2) | ((0x2000 - d) & 0x1f00) >> 5);
      } else {
        ret.write<d8, d8>(((0x2000 - d) & 0x1f00) >> 5, cmd.len - 1);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh, d16, d8>({2, 1}, 0xf000, 0);
  ret.trim();
  assert(adr == input.size());
  assert(dp.optimal_cost() + 2 + 3 * 8 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
