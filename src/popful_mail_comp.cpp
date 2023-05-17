#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> popful_mail_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum tag { uncomp, uncompl, rle, rlel, lz };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x1003);
  auto c1 = dp.c<1>(0x1fff);
  auto c1_31 = dp.c<1, 31>(input.size());

  size_t rlen = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 0x1f, c1, 1, uncomp);
    dp.update(i, 2, 0x1fff, c1, 2, uncompl);
    rlen = encode::run_length_r(input, i, rlen);
    dp.update(i, 4, 0x13, rlen, c0, 2, rle);
    dp.update(i, 4, 0x1003, rlen, c0, 3, rlel);
    const auto res_lz = lz_helper.find(i, 0x1fff, 4);
    for (size_t r = 4; r < std::min<size_t>(res_lz.len + 1, 4 + 31); ++r) {
      dp.update(i, r, input.size(), res_lz, c1_31, 2 + (r % 31 + 23) / 31, lz);
    }
    c0.update(i); c1_31.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
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
  assert(dp.optimal_cost() + 3 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
