#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> mujintou_monogatari_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  enum tag { uncomp, lz, lzl };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(3 + 0x0f + 0xff);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    const auto res_lz = lz_helper.find(i, 0x0fff, 3);
    dp.update(i, 3, 3 + 0x0e, res_lz, c0, 17, lz);
    dp.update(i, 3 + 0x0f, 3 + 0x0f + 0xff, res_lz, c0, 25, lzl);
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(3);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: ret.write<b1, d8>(false, input[adr]); break;
    case lz: ret.write<b1, d16>(true, d << 4 | (cmd.len - 3)); break;
    case lzl: ret.write<b1, d16, d8>(true, d << 4 | 0x000f, cmd.len - 3 - 0x0f); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 3 * 8 == ret.bit_length());
  ret[0] = 0x0c;
  write16(ret.out, 1, input.size());
  return ret.out;
}

} // namespace sfc_comp
