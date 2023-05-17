#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> wizardry5_comp_1(std::span<const uint8_t> input) {
  enum tag { uncomp0, uncomp, lz };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());

  static constexpr auto lz_lens = std::to_array<size_t>({
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 32, 64, 128, 256
  });
  static constexpr auto len_code = inverse_map<lz_lens.back() + 1>(lz_lens);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    if (input[i] == 0) dp.update(i, 1, 5, uncomp0);
    else dp.update(i, 1, 9, uncomp);
    dp.update(i, lz_lens, lz_helper.find(i, 0x200, lz_lens.front()), constant<14>(), lz);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp0: ret.write<bnh>({5, 0x10}); break;
    case uncomp: ret.write<bnh>({9, input[adr]}); break;
    case lz: ret.write<bnh>({14, (1 << 13) | (len_code[cmd.len] + 1) << 9 | (d - 1)}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({9, 0});
  assert(dp.optimal_cost() + 9 == ret.bit_length());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
