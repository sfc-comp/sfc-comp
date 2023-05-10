#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> legend_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum tag { uncomp, lz2, lzs, lzm, lzl };
  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find(i, 0x1000, 2);
    dp.update_lz(i, 2, 2, res_lz, Constant<14>(), lz2);
    dp.update_lz(i, 3, 10, res_lz, Constant<18>(), lzs);
    dp.update_lz(i, 11, 18, res_lz, Constant<19>(), lzm);
    dp.update_lz(i, 19, 274, res_lz, Constant<24>(), lzl);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(8);
  writer_b8_h flags;

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      flags.write<b1>(false);
      ret.write<d8>(input[adr]);
    } break;
    case lz2:
    case lzs:
    case lzm:
    case lzl: {
      flags.write<b1>(true);
      if (cmd.type == lz2) {
        flags.write<bnh>({1, 1});
      } else if (cmd.type == lzs) {
        flags.write<bnh>({5, 8 | (cmd.len - 3)});
      } else if (cmd.type == lzm) {
        flags.write<bnh>({6, 8 | (cmd.len - 11)});
      } else {
        flags.write<bnh>({3, 0});
        ret.write<d8>(cmd.len - 19);
      }
      const size_t d = adr - cmd.lz_ofs;
      ret.write<d8>((d - 1) >> 4);
      flags.write<bnh>({4, (d - 1) & 0x0f});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 8 * 8 == ret.size() * 8 + flags.bit_length());
  write32b(ret.out, 0, input.size());
  write32b(ret.out, 4, ret.size() - 8);
  ret.extend(flags);
  return ret.out;
}

} // namespace sfc_comp
