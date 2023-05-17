#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> legend_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum method { uncomp, lz };
  using tag = tag_l<method>;

  static constexpr auto lens = to_vranges({
    {0x002,  1, 0b1},
    {0x003,  5, 0b01'000},
    {0x00b,  6, 0b001'000},
    {0x013, 11, 0b000'00000000}
  }, 0x100); // cf. $82:8868

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(lens.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, {uncomp, 0});
    const auto res_lz = lz_helper.find(i, 0x1000, 2);
    dp.update(i, lens, res_lz, c0, 13, [&](size_t li) -> tag { return {lz, li}; });
    c0.update(i);
  }

  using namespace data_type;
  writer ret(8);
  writer_b8_h flags;

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      flags.write<b1>(false);
      ret.write<d8>(input[adr]);
    } break;
    case lz: {
      const auto& l = lens[li];
      const size_t v = l.val + (cmd.len - l.min);
      flags.write<b1>(true);
      if (l.bitlen < 8) {
        flags.write<bnh>({l.bitlen, v});
      } else {
        flags.write<bnh>({l.bitlen - 8, v >> 8});
        ret.write<d8>(v & 0xff);
      }
      const size_t d = adr - cmd.lz_ofs();
      ret.write<d8>((d - 1) >> 4);
      flags.write<bnh>({4, (d - 1) & 0x0f});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 8 * 8 == ret.size() * 8 + flags.bit_length());
  write32b(ret.out, 0, input.size());
  write32b(ret.out, 4, ret.size() - 8);
  ret.extend(flags);
  return ret.out;
}

} // namespace sfc_comp
