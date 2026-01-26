#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> asuka_kenzan_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x800000);

  enum tag { uncomp0, uncomp1, lzs, lzl };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x103);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    if (input[i] & 0x80) dp.update(i, 1, 9, uncomp1);
    else dp.update(i, 1, 8, uncomp0);
    dp.update(i, 3, 0x12, lz_helper.find(i, 0x800, 3), c0, 18, lzs);
    dp.update(i, 4, 0x103, lz_helper.find(i, 0x8000, 3), c0, 26, lzl);
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_pre_l ret(5); ret.write<none>(none());
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp0: ret.write<d8>(input[adr]); break;
    case uncomp1: ret.write<d8, b1>(input[adr], false); break;
    case lzs: ret.write<d8, b1, b1, d8>((0x1000 - d) >> 4, true, false, ((0x1000 - d) & 0x0f) << 4 | (cmd.len - 3)); break;
    case lzl: ret.write<d8, b1, b1, d8, d8>((0x10000 - d) >> 8, true, true, (0x10000 - d) & 0xff, (cmd.len - 4)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret[0] = 0x34;
  write32(ret.out, 1, input.size());
  assert(adr == input.size());
  assert(dp.optimal_cost() + 5 * 8 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
