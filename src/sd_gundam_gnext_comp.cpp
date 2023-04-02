#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sd_gundam_gnext_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);
  enum CompType {
    uncomp, lz
  };
  const size_t pad = 0x22;
  std::vector<uint8_t> input2(input.size() + pad);
  std::copy(input.begin(), input.end(), input2.begin() + pad);

  lz_helper lz_helper(input2);
  sssp_solver<CompType> dp(input2.size());

  for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);
  for (size_t i = 0; i < pad; ++i) dp[i + 1].cost = 0;

  for (size_t i = pad; i < input2.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x1000);
    dp.update_k<2>(i, 3, 0x21, res_lz.len, Constant<17>(), lz, res_lz.ofs);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  size_t adr = pad;
  ret.write<d16>(0);
  for (const auto cmd : dp.commands(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<b1l, d8>(true, input2[adr]); break;
    case lz: {
      assert(cmd.len & 1);
      uint16_t d = (cmd.lz_ofs - pad) & 0x0fff;
      ret.write<b1l, d8, d8>(false, d & 0x00ff, (d & 0x0f00) >> 8 | (cmd.len - 3) << 3);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert((dp.total_cost() + 7) / 8 + 2 == ret.size());
  assert(adr - pad == input.size());

  if (ret.size() < input.size() + 2) {
    write16(ret.out, 0, (input.size() - 1) | 0x8000);
  } else {
    write16(ret.out, 0, (input.size() - 1) | 0x0000);
    ret.out.resize(input.size() + 2);
    std::copy(input.begin(), input.end(), ret.out.begin() + 2);
  }
  return ret.out;
}

} // namespace sfc_comp
