#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> power_piggs_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);

  enum tag { uncomp, lz };
  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.c<0>(0x12);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    dp.update(i, 3, 18, lz_helper.find(i, 0x8000, 3), c0, 20, lz);
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    switch (cmd.type) {
    case uncomp: {
      ret.write<b1, bnh>(true, {8, input[adr]});
    } break;
    case lz: {
      ret.write<b1>(false);
      ret.write<bnh, bnh>({15, cmd.lz_ofs() + 1}, {4, cmd.len - 3});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({16, 0x0000});
  assert(adr == input.size());
  assert(dp.optimal_cost() + 16 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
