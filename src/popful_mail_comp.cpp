#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> popful_mail_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum CompType {
    uncomp, uncompl, rle, rlel, lz
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 0x1f, Linear<1, 1>(), uncomp);
    dp.update(i, 2, 0x1fff, Linear<1, 2>(), uncompl);
    rlen = encode::run_length(input, i, rlen);
    dp.update(i, 4, 0x13, rlen, Constant<2>(), rle);
    dp.update(i, 4, 0x1003, rlen, Constant<3>(), rlel);
    auto res_lz = lz_helper.find_best(i, 0x1fff);
    for (size_t j = 0; j < 31; ++j) {
      if (j > res_lz.len) break;
      size_t beg = (j < 4) ? j + 31 : j;
      dp.update_k<31>(i, beg, 0x8000, res_lz.len, LinearQ<1, 93 - 8, 31>(), lz, res_lz.ofs);
    }
    lz_helper.add_element(i);
  }
  using namespace data_type;
  writer ret; ret.write<d16b>(0);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(cmd.len, {cmd.len, &input[adr]}); break;
    case uncompl: ret.write<d16b, d8n>(0x2000 + cmd.len, {cmd.len, &input[adr]}); break;
    case rle: ret.write<d8, d8>(0x40 | (cmd.len - 4), input[adr]); break;
    case rlel: ret.write<d16b, d8>(0x5000 | (cmd.len - 4), input[adr]); break;
    case lz: {
      if (cmd.len >= 8) {
        ret.write<d16b>(0xe000 | d);
        size_t l = cmd.len - 7;
        for (; l >= 0x1f; l -= 0x1f) ret.write<d8>(0x7f);
        if (l > 0) ret.write<d8>(0x60 + l);
      } else {
        ret.write<d16b>(0x8000 | (cmd.len - 4) << 13 | d);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d8>(0);
  write16(ret.out, 0, ret.out.size() - 2);
  assert(dp.total_cost() + 3 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
