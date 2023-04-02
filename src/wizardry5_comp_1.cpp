#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> wizardry5_comp_1(std::span<const uint8_t> input) {
  enum CompType {
    uncomp0, uncomp, lz
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  static constexpr std::array<size_t, 15> lz_len_table = {
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 32, 64, 128, 256
  };
  std::vector<uint8_t> code(257);
  for (size_t i = 0; i < lz_len_table.size(); ++i) code[lz_len_table[i]] = i + 1;

  for (size_t i = 0; i < input.size(); ++i) {
    if (input[i] == 0) dp.update(i, 1, 1, Constant<5>(), uncomp0);
    else dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x200);
    dp.update_lz_table(i, lz_len_table, res_lz, Constant<14>(), lz);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp0: ret.write<b8hn_h>({5, 0x10}); break;
    case uncomp: ret.write<b8hn_h>({9, input[adr]}); break;
    case lz: ret.write<b8hn_h>({14, (1 << 13) | code[cmd.len] << 9 | (d - 1)}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b8hn_h>({9, 0});
  assert((dp.total_cost() + 7 + 9) / 8 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
