#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> legend_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum CompType {
    uncomp, lz2, lzs, lzm, lzl
  };
  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x1000);
    dp.update_lz(i, 2, 2, res_lz, Constant<14>(), lz2);
    dp.update_lz(i, 3, 10, res_lz, Constant<18>(), lzs);
    dp.update_lz(i, 11, 18, res_lz, Constant<19>(), lzm);
    dp.update_lz(i, 19, 274, res_lz, Constant<24>(), lzl);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret; ret.write<d32b, d32b>(0, 0);
  writer_b flags;

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      flags.write<b1h>(false);
      ret.write<d8>(input[adr]);
    } break;
    case lz2:
    case lzs:
    case lzm:
    case lzl: {
      flags.write<b1h>(true);
      if (cmd.type == lz2) {
        flags.write<b8hn_h>({1, 1});
      } else if (cmd.type == lzs) {
        flags.write<b8hn_h>({5, 8 | (cmd.len - 3)});
      } else if (cmd.type == lzm) {
        flags.write<b8hn_h>({6, 8 | (cmd.len - 11)});
      } else {
        flags.write<b8hn_h>({3, 0});
        ret.write<d8>(cmd.len - 19);
      }
      const size_t d = adr - cmd.lz_ofs - 1;
      ret.write<d8>(d >> 4);
      flags.write<b8hn_h>({4, d & 0x0f});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write32b(ret.out, 0, input.size());
  write32b(ret.out, 4, ret.size() - 8);
  std::copy(flags.out.begin(), flags.out.end(), std::back_inserter(ret.out));
  assert(adr == input.size());
  assert((dp.total_cost() + 7) / 8 + 8 == ret.size());
  return ret.out;
}

} // namespace sfc_comp
