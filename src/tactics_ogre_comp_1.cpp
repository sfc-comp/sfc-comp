#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> tactics_ogre_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum tag { uncomp, rle0, lzs, lzl };

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 64, Linear<1, 1>(), uncomp);
    rlen = encode::run_length(input, i, rlen);
    if (input[i] == 0) dp.update(i, 2, 0x21, rlen, Constant<1>(), rle0);
    auto res_lzs = lz_helper.find(i, 0x0800, 3);
    dp.update_lz(i, 3, 0x12, res_lzs, Constant<2>(), lzs);
    auto res_lzl = lz_helper.find(i, 0x4000, 4);
    dp.update_lz(i, 4, 0x43, res_lzl, Constant<3>(), lzl);
    lz_helper.add_element(i);
  }
  using namespace data_type;
  writer ret(2);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x40 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case rle0: ret.write<d8>(0x20 + cmd.len - 2); break;
    case lzs: ret.write<d16b>(0x8000 | ((cmd.len - 3) << 11) | (d - 1)); break;
    case lzl: ret.write<d24b>(
      0x100000 | ((cmd.len - 4) & 0x0f) << 16
               | ((cmd.len - 4) & 0x30) << 14 | (d - 1)); break;
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
