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
  writer_b8_h ret(2);
  size_t adr = pad;
  for (const auto cmd : dp.commands(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<bnh>({9, size_t(0x100 | input2[adr])}); break;
    case lz: ret.write<bnh, bnh>({9, (cmd.lz_ofs - pad - 0x11) & 0xff}, {4, cmd.len - 2}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());

  assert(adr - pad == input.size());
  assert(dp.total_cost() + 2 * 8 == ret.bit_length());

  return ret.out;
}

std::vector<uint8_t> sotsugyou_bangai_hen_comp(std::span<const uint8_t> input) {
  auto ret = slap_stick_comp(input);
  for (size_t i = 0; i < ret.size() - 2; ++i) ret[i] = ret[i + 2];
  ret.resize(ret.size() - 2);
  return ret;
}

} // namespace sfc_comp
