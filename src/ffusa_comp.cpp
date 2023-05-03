#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> ffusa_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum tag { uncomp, lz };

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());
  using cost_type = typename sssp_solver<tag>::cost_type;

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update<std::greater<cost_type>>(i, 1, 0x0f, Linear<1, 1>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x0100);
    if (i > 0) {
      if (dp[i].type == uncomp) {
        dp.update_lz(i, 3, 0x11, res_lz, Constant<1>(), lz);
      } else {
        dp.update_lz(i, 3, 0x11, res_lz, Constant<2>(), lz);
      }
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer raw, ret; ret.write<d16>(0);
  size_t adr = 0;

  bool uncomped = false;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      ret.write<d8>(cmd.len);
      raw.write<d8n>({cmd.len, &input[adr]});
      uncomped = true;
    } break;
    case lz: {
      if (uncomped) {
        ret.out.back() |= (cmd.len - 2) << 4;
      } else {
        ret.write<d8>((cmd.len - 2) << 4);
      }
      ret.write<d8>((adr - cmd.lz_ofs) - 1);
      uncomped = false;
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d8>(0);
  write16(ret.out, 0, ret.size() - 2);
  ret.extend(raw);

  assert(dp.total_cost() + 3 == ret.size());
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
