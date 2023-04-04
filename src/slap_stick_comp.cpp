#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> slap_stick_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum CompType {
    uncomp, lz
  };

  static constexpr size_t pad = 0x11;

  std::vector<uint8_t> input2(input.size() + pad);
  std::copy(input.begin(), input.end(), input2.begin() + pad);
  for (size_t i = 0; i < pad; ++i) input2[i] = 0x20;

  lz_helper lz_helper(input2);
  sssp_solver<CompType> dp(input2.size());

  for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);
  for (size_t i = 0; i < pad; ++i) dp[i + 1].cost = 0;

  for (size_t i = pad; i < input2.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x100);
    dp.update_lz(i, 2, 0x11, res_lz, Constant<13>(), lz);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret; ret.write<d16>(0);
  size_t adr = pad;
  for (const auto cmd : dp.commands(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<b8hn_h>({9, size_t(0x100 | input2[adr])}); break;
    case lz: ret.write<b8hn_h, b8hn_h>({9, (cmd.lz_ofs - pad * 2) & 0xff}, {4, cmd.len - 2}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());

  assert((dp.total_cost() + 7) / 8 + 2 == ret.size());
  assert(adr - pad == input.size());

  return ret.out;
}

} // namespace sfc_comp
