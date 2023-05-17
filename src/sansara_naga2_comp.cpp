#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sansara_naga2_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  enum tag { uncomp, lzs, lzl, lzll };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x4000);
  auto c1 = dp.c<1>(0x3f);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 0x3f, c1, 1, uncomp);
    dp.update(i, 2, 5, lz_helper.find(i, 0x10, 2), c0, 1, lzs);
    const auto res_lzl = lz_helper.find(i, 0x400, 2);
    dp.update(i, 1, 0x20, res_lzl, c0, 2, lzl);
    dp.update(i, 1, 0x4000, res_lzl, c0, 4, lzll);
    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret(6);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case lzll: ret.write<d8, d24b>(0x00,  (cmd.len - 1) << 10 | (d - 1)); break;
    case uncomp: ret.write<d8, d8n>(cmd.len, {cmd.len, &input[adr]}); break;
    case lzs: ret.write<d8>(0x40 | ((cmd.len - 2) << 4) | (d - 1)); break;
    case lzl: ret.write<d16b>(0x8000 | ((cmd.len - 1) << 10) | (d - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret[0] = 'P'; ret[1] = '2';
  write16(ret.out, 2, input.size());
  write16(ret.out, 4, ret.out.size() - 6);
  assert(dp.optimal_cost() + 6 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
