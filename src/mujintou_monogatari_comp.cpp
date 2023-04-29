#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> mujintou_monogatari_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  enum CompType { uncomp, lz, lzl };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x0fff);
    dp.update_lz(i, 3, 3 + 0x0e, res_lz, Constant<17>(), lz);
    dp.update_lz(i, 3 + 0x0f, 3 + 0x0f + 0xff, res_lz, Constant<25>(), lzl);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret(2);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<b1, d8>(false, input[adr]); break;
    case lz: ret.write<b1, d16>(true, d << 4 | (cmd.len - 3)); break;
    case lzl: ret.write<b1, d16, d8>(true, d << 4 | 0x000f, cmd.len - 3 - 0x0f); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.total_cost() + 2 * 8 == ret.bit_length());
  write16(ret.out, 0, input.size());
  return ret.out;
}

} // namespace sfc_comp
