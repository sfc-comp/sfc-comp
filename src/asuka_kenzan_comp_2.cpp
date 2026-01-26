#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> asuka_kenzan_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz, lz2 };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  8, 0b00'000000},
    {0x0041, 10, 0b0100'000000},
    {0x00c1, 11, 0b01100'000000},
    {0x01c1, 12, 0b100000'000000},
    {0x0581, 13, 0b1011110'000000},
    {0x0c01, 14, 0b11110000'000000}
  }, 0x1000);

  static constexpr auto ofs2_tab = create_array<vrange, ofs_tab.size()>([](size_t i) {
    const auto o = ofs_tab[i];
    return vrange(((o.min - 1) >> 4) + 1, o.max >> 4, o.bitlen - 4, o.val >> 4);
  });

  static constexpr auto len_tab = to_vranges({
              // 0b010 (len == 2)
    {0x0003,  2, 0b00},
    {0x0004,  3, 0b011},
    {0x0006,  4, 0b1010},
    {0x0009,  5, 0b11010},
    {0x000a,  6, 0b11011'0},
    {0x000c,  7, 0b11100'00},
    {0x0010,  8, 0b11101'000},
    {0x0018, 10, 0b111100'0000},
    {0x0028, 11, 0b111101'00000},
    {0x0048, 12, 0b111110'000000},
    {0x0088, 14, 0b1111110'0000000},
    {0x0108, 15, 0b1111111'00000000},
  }, 0x0206);

  static constexpr auto len2_tab = to_vranges({
    {0x0002,  3, 0b010}
  }, 0x0002);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.c<0>(len_tab.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, {uncomp, 0, 0});
    dp.update_matrix(i, ofs2_tab, len2_tab, c0, 1,
      [&](size_t oi) { return lz_helper.find(i, ofs2_tab[oi].max, len2_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz2, oi, li}; }
    );
    dp.update_matrix(i, ofs_tab, len_tab, c0, 1,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_l ret; ret.write<bnl>({16, 0x0600});

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    const size_t d = adr - cmd.lz_ofs();
    switch (tag) {
    case uncomp: {
      ret.write<b1, bnl>(false, {8, input[adr]});
    } break;
    case lz2: {
      const auto& l = len2_tab[li];
      const auto& o = ofs2_tab[oi];
      ret.write<b1>(true);
      ret.write<bnh>({l.bitlen, ~l.val});
      ret.write<bnh>({o.bitlen - 2, ~((o.val + (d - o.min)) >> 2)});
      ret.write<bnl>({2, d - o.min});
    } break;
    case lz: {
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[oi];
      const size_t s = (li < 4) ? 0 : (li - 3);
      ret.write<b1>(true);
      ret.write<bnh>({l.bitlen - s, ~((l.val + (cmd.len - l.min)) >> s)});
      ret.write<bnl>({s, cmd.len - l.min});
      ret.write<bnh>({o.bitlen - 6, ~((o.val + (d - o.min)) >> 6)});
      ret.write<bnl>({6, d - o.min});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1>(true);
  ret.write<bnh, bnl>({7, 0}, {8, 0xff});
  assert(adr == input.size());
  assert(dp.optimal_cost() + 16 + 1 + 15 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
