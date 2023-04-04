#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> seiken_densetsu_2_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum CompType {
    uncomp, lz
  };

  std::vector<uint8_t> best;
  for (size_t ty = 0; ty < 6; ++ty) {
    const size_t min_len = 3;
    const size_t max_len = min_len + (0x007f >> (5 - ty));
    const size_t max_ofs = (0x2000 >> ty);

    lz_helper lz_helper(input);
    sssp_solver<CompType> dp(input.size());

    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 0x80, Linear<1, 1>(), uncomp);
      auto res_lz = lz_helper.find_best(i, max_ofs);
      dp.update_lz(i, min_len, max_len, res_lz, Constant<2>(), lz);
      lz_helper.add_element(i);
    }

    using namespace data_type;
    writer ret; ret.write<d16, d16>(ty, 0);
    size_t adr = 0;

    for (const auto cmd : dp.commands()) {
      switch (cmd.type) {
      case uncomp: {
        ret.write<d8, d8n>(cmd.len - 1, {cmd.len, &input[adr]});
      } break;
      case lz: {
        ret.write<d16b>((adr - cmd.lz_ofs - 1) | (cmd.len - 3) << (13 - ty) | 0x8000);
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    write16b(ret.out, 1, input.size());
    assert(dp.total_cost() + 4 == ret.size());
    if (best.size() == 0 || ret.size() < best.size()) {
      best = std::move(ret.out);
    }
  }
  return best;
}

} // namespace sfc_comp
