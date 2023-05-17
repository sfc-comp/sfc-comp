#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> rayearth_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ulens = to_vranges({
    {0x01,  8, 0x00},
    {0x40, 16, 0x3f00}
  }, 0xff);

  static constexpr auto lens = to_vranges({
    {0x03,  6, 0x00},
    {0x42, 14, 0x3f00}
  }, 0xff);

  static constexpr auto offsets = to_vranges({
    {0x0001, 10, 0x4000 + 1, 0xc0ff},
    {0x02ff, 18, 0xc0ff00,   0x0000ff}
  }, 0x03fe);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(lens.back().max);
  auto c8 = dp.c<8>(ulens.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, ulens, c8, 0, [&](size_t li) -> tag { return {uncomp, 0, li}; });
    dp.update_matrix(i, offsets, lens, c0, 0,
      [&](size_t oi) { return lz_helper.find(i, offsets[oi].max, lens.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer ret(3);

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    const size_t d = adr - cmd.lz_ofs();
    switch (tag) {
    case uncomp: {
      const auto& u = ulens[li];
      const size_t v = u.val + (cmd.len - u.min);
      if (u.bitlen <= 8) ret.write<d8>(v);
      else ret.write<d16b>(v);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const auto& l = lens[li];
      const auto& o = offsets[oi];
      const size_t lv = l.val + (cmd.len - l.min);
      const size_t ov = masked_add(o.val, d - o.min, o.mask);
      ret.write<d8>(lv >> (l.bitlen - 6) | ov >> (o.bitlen - 2));
      if (l.bitlen > 6) ret.write<d8>(lv & 0xff);
      if (o.bitlen == 10) ret.write<d8>(ov & 0xff);
      else ret.write<d16b>(ov & 0xffff);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret[0] = 0;
  write16(ret.out, 1, input.size());
  assert(dp.optimal_cost() + 3 * 8 == ret.size() * 8);
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
