#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sailor_moon_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum tag { uncomp, lzs, lzls, lzll };
  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lzs = lz_helper.find_best(i, 0x100);
    dp.update_lz(i, 2, 5, res_lzs, Constant<12>(), lzs);
    auto res_lzl = lz_helper.find_best(i, 0x2000);
    dp.update_lz(i, 3, 9, res_lzl, Constant<18>(), lzls);
    dp.update_lz(i, 10, 256, res_lzl, Constant<26>(), lzll);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16_pre_l ret;
  ret.write<none>(none());

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
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
  assert(dp.total_cost() + 2 + 3 * 8 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
