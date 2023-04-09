#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> danzarb_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum CompType {
    uncomp, lz
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 8, Linear<8, 4>(), uncomp);
    if (i > 0) {
      auto res_lz = lz_helper.find_best(i, input.size());
      size_t c = ilog2(2 * i - 1);
      dp.update_lz(i, 1, 32, res_lz, Constant<1 + 5>(), lz, dp[i].cost + c);
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret; ret.write<d16>(0);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      ret.write<b1h, b8hn_h>(false, {3, cmd.len & 7});
      for (size_t i = 0; i < cmd.len; ++i) ret.write<b8hn_h>({8, input[adr + i]});
    } break;
    case lz: {
      ret.write<b1h, b8hn_h, b8hn_h>(true, {ilog2(2 * adr - 1), cmd.lz_ofs}, {5, cmd.len & 0x1f});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert((dp.total_cost() + 7) / 8 + 2 == ret.size());
  assert(adr == input.size());
  write16(ret.out, 0, input.size());
  return ret.out;
}

} // namespace sfc_comp
