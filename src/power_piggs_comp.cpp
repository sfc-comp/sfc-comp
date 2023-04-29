#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> power_piggs_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);

  enum CompType {
    uncomp, lz
  };
  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x8000);
    dp.update_lz(i, 3, 18, res_lz, Constant<20>(), lz);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      ret.write<b1, bnh>(true, {8, input[adr]});
    } break;
    case lz: {
      ret.write<b1>(false);
      ret.write<bnh, bnh>({15, cmd.lz_ofs + 1}, {4, cmd.len - 3});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({16, 0x0000});
  assert(adr == input.size());
  assert(dp.total_cost() + 16 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
