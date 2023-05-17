#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> tactics_ogre_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum tag { uncomp, rle0, lzs, lzl };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x43);
  auto c1 = dp.c<1>(0x40);

  size_t rlen = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 64, c1, 1, uncomp);
    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] == 0) dp.update(i, 2, 0x21, rlen, c0, 1, rle0);
    dp.update(i, 3, 0x12, lz_helper.find(i, 0x0800, 3), c0, 2, lzs);
    dp.update(i, 4, 0x43, lz_helper.find(i, 0x4000, 4), c0, 3, lzl);
    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case lzl: ret.write<d24b>(0x100000 | ((cmd.len - 4) & 0x0f) << 16
                                       | ((cmd.len - 4) & 0x30) << 10 | (d - 1)); break;
    case rle0: ret.write<d8>(0x20 + cmd.len - 2); break;
    case uncomp: ret.write<d8, d8n>(0x40 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case lzs: ret.write<d16b>(0x8000 | ((cmd.len - 3) << 11) | (d - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  assert(dp.optimal_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
