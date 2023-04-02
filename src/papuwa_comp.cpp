#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> papuwa_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum CompType {
    uncomp, uncompl,
    lzs, lzm, lzls, lzlm, lzll
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 16, Linear<1, 1>(), uncomp);
    dp.update(i, 17, 0x1010, Linear<1, 2>(), uncompl);
    auto res_lzs = lz_helper.find_best(i, 0x10);
    dp.update_lz(i, 3, 6, res_lzs, Constant<1>(), lzs);
    auto res_lzm = lz_helper.find_best(i, 0x400);
    dp.update_lz(i, 7, 22, res_lzm, Constant<2>(), lzm);
    auto res_lzl = lz_helper.find_best(i, 0x1000);
    dp.update_lz(i, 3, 6, res_lzl, Constant<2>(), lzls);
    dp.update_lz(i, 7, 0x206, res_lzl, Constant<3>(), lzlm);
    dp.update_lz(i, 0x207, 0x8206, res_lzl, Constant<4>(), lzll);
    lz_helper.add_element(i);
  }
  using namespace data_type;
  writer ret; ret.write<d16b>(0);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0xe0 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case uncompl: ret.write<d16b, d8n>(0xf800 + cmd.len - 17, {cmd.len, &input[adr]}); break;
    case lzs: ret.write<d8>((cmd.len - 3) | (d - 1) << 2); break;
    case lzm: ret.write<d16b>(0x8000 | (cmd.len - 7) << 10 | (d - 1)); break;
    case lzls: ret.write<d16b>(
      0x4000 | ((d - 1) & 0x0f00) << 2 | (cmd.len - 3) << 8
             | ((d - 1) & 0x00ff)); break;
    case lzlm: ret.write<d24b>(
      0xc00000 | ((cmd.len - 7) & 0x100) << 12
               | ((cmd.len - 7) & 0x00f) << 16
               | ((cmd.len - 7) & 0x0f0) << 8 | (d - 1)); break;
    case lzll: ret.write<d32b>(
      0xf0000000 | ((cmd.len - 519) & 0x7000) << 12
                 | ((cmd.len - 519) & 0x00ff) << 16
                 | ((cmd.len - 519) & 0x0f00) << 4
                 | (d - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  assert(dp.total_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
