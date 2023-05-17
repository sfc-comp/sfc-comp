#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> ffusa_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum tag { uncomp, lz, lz_none };

  lz_helper lz_helper(input, true);
  solver<tag> dp0(input.size()), dp1(input.size());
  auto c0_0 = dp0.c<0>(0x11);
  auto c1_1 = dp1.c<1>(0x0f);

  for (size_t i = input.size(); ; ) {
    dp0.update_c(i, 0, dp1[i].cost + 1, uncomp);
    c0_0.update(i);
    dp1.update_c(i, 0, dp0[i].cost, lz_none);
    c1_1.update(i);
    if (i-- == 0) break;
    lz_helper.reset(i);
    dp0.update(i, 1, 0x0f, c1_1, 1, uncomp);
    dp1.update(i, 3, 0x11, lz_helper.find(i, 0x0100, 3), c0_0, 1, lz);
  }

  using namespace data_type;
  writer ret(2); writer raw;

  size_t adr = 0;
  for (size_t curr = 0; adr < input.size(); curr ^= 1) {
    const auto& cmd = (curr == 0) ? dp0[adr] : dp1[adr];
    switch (cmd.type) {
    case uncomp: {
      ret.write<d8>(cmd.len);
      raw.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz_none: break;
    case lz: {
      ret[ret.size() - 1] |= (cmd.len - 2) << 4;
      ret.write<d8>((adr - cmd.lz_ofs()) - 1);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d8>(0);
  write16(ret.out, 0, ret.size() - 2);
  ret.extend(raw);

  assert(dp0.optimal_cost() + 3 == ret.size());
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
