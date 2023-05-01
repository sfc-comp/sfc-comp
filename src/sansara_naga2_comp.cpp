#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sansara_naga2_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  enum tag { uncomp, lzs, lzl, lzll };

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 0x3f, Linear<1, 1>(), uncomp);
    auto res_lzs = lz_helper.find_best(i, 0x10);
    dp.update_lz(i, 2, 5, res_lzs, Constant<1>(), lzs);
    auto res_lzl = lz_helper.find_best(i, 0x400);
    dp.update_lz(i, 1, 0x20, res_lzl, Constant<2>(), lzl);
    dp.update_lz(i, 1, 0x4000, res_lzl, Constant<4>(), lzll);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(6);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(cmd.len, {cmd.len, &input[adr]}); break;
    case lzs: ret.write<d8>(0x40 | ((cmd.len - 2) << 4) | (d - 1)); break;
    case lzl: ret.write<d16b>(0x8000 | ((cmd.len - 1) << 10) | (d - 1)); break;
    case lzll: ret.write<d8, d24b>(0x00,  (cmd.len - 1) << 10 | (d - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret[0] = 'P';
  ret[1] = '2';
  write16(ret.out, 2, input.size());
  write16(ret.out, 4, ret.out.size() - 6);
  assert(dp.total_cost() + 6 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
