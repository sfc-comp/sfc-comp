#include "algorithm.hpp"
#include "encode.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> front_mission_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum CompType {
    uncomp, rle, lz
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  encode::lz_data res_lz_curr = {}, res_lz_next = {};
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz1 = lz_helper.find_best(i, 1);
    dp.update(i, 3, 0x82, res_lz1.len, Constant<9>(), rle);
    res_lz_curr = res_lz_next;
    if (i + 1 < input.size()) {
      res_lz_next = lz_helper.find_best(i + 1, 0xfff);
    }
    dp.update_k<2>(
      i, 4, 0x12, res_lz_curr.len, Constant<17>(), lz, res_lz_curr.ofs);
    lz_helper.add_element(i);
  }
  using namespace data_type;
  writer_b ret; ret.write<d16>(0);
  size_t adr = 0, num_codes = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<b1l, d8>(false, input[adr]); break;
    case rle: ret.write<b1l, d8>(true, (cmd.len - 3) << 1); break;
    case lz: ret.write<b1l, d16b>(true,
      (d & 0xf00) << 4 | (cmd.len - 4) << 8 | 0x100 | (d & 0x0ff)); break;
    default: assert(0);
    }
    adr += cmd.len;
    num_codes += 1;
  }
  write16(ret.out, 0, num_codes);
  assert((dp.total_cost() + 7) / 8 + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
