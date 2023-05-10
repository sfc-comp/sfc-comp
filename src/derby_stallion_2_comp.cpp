#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> derby_stallion_2_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  9, 0b0'00000000 + 1},
    {0x0100, 12, 0b1'00000000000},
  }, 0x08ff);

  static constexpr auto len_tab = to_vranges({
    {0x0002,  1, 0b1},
    {0x0003,  2, 0b01},
    {0x0004,  5, 0b001'00},
    {0x0008,  7, 0b0001'000},
    {0x0010,  9, 0b0000'00000},
  }, 0x002f);

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      1
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(2); writer_b8_l flags;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      flags.write<b1>(true);
      ret.write<d8>(input[adr]);
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const auto& o = ofs_tab[cmd.type.oi];
      const auto& l = len_tab[cmd.type.li];
      size_t v = (d - o.min) + o.val;
      flags.write<b1>(false);
      flags.write<bnh>({o.bitlen - 8, v >> 8});
      ret.write<d8>(v & 0xff);
      flags.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  flags.write<b1, b1>(false, false);
  ret.write<d8>(0);
  assert(adr == input.size());
  assert(dp.optimal_cost() + 2 * 8 + 10 == flags.bit_length() + ret.size() * 8);
  write16(ret.out, 0, ret.size() - 2);
  ret.extend(flags);
  return ret.out;
}

} // namespace sfc_comp
