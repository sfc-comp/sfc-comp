#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> papuwa_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum tag {
    uncomp, uncompl,
    lzs, lzm, lzls, lzlm, lzll
  };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x8206);
  auto c1 = dp.c<1>(0x11 + 0x03ff);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 0x10, c1, 1, uncomp);
    dp.update(i, 0x11, 0x11 + 0x03ff, c1, 2, uncompl);
    dp.update(i, 3, 6, lz_helper.find(i, 0x10, 3), c0, 1, lzs);
    dp.update(i, 7, 22, lz_helper.find(i, 0x400, 7), c0, 2, lzm);
    const auto res_lzl = lz_helper.find(i, 0x1000, 3);
    dp.update(i, 3, 6, res_lzl, c0, 2, lzls);
    dp.update(i, 7, 0x206, res_lzl, c0, 3, lzlm);
    dp.update(i, 0x207, 0x8206, res_lzl, c0, 4, lzll);
    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case lzs: ret.write<d8>((cmd.len - 3) | (d - 1) << 2); break;
    case lzls: ret.write<d16b>(0x4000 | ((d - 1) & 0x0f00) << 2 | (cmd.len - 3) << 8
                                      | ((d - 1) & 0x00ff)); break;
    case lzm: ret.write<d16b>(0x8000 | (cmd.len - 7) << 10 | (d - 1)); break;
    case lzlm: ret.write<d24b>(0xc00000 | ((cmd.len - 7) & 0x100) << 12
                                        | ((cmd.len - 7) & 0x00f) << 16
                                        | ((cmd.len - 7) & 0x0f0) << 8 | (d - 1)); break;
    case uncomp: ret.write<d8, d8n>(0xe0 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case lzll: ret.write<d32b>(0xf0000000 | ((cmd.len - 519) & 0x7000) << 12
                                          | ((cmd.len - 519) & 0x00ff) << 16
                                          | ((cmd.len - 519) & 0x0f00) << 4
                                          | (d - 1)); break;
    case uncompl: ret.write<d16b, d8n>(0xf800 + cmd.len - 0x11, {cmd.len, &input[adr]}); break;
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
