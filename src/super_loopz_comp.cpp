#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_loopz_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0xffff);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  7, 0b10'00000 + 1},
    {0x0020, 10, 0b0'000000000},
    {0x0220, 16, 0b11'00000000000000},
  }, 0x421f);
  static constexpr std::array<size_t, ofs_tab.size()> ofs_prefix = {2, 1, 2};

  static constexpr auto len_tab = to_vranges({
    {0x0002,  2, 0b0'0},
    {0x0004,  4, 0b10'00},
    {0x0008,  7, 0b110'0000},
    {0x0017, 11, 0b111'00000000},
  }, 0x0116);
  static constexpr std::array<size_t, len_tab.size()> len_prefix = {1, 2, 3, 3};

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(len_tab.back().max);
  auto c8 = dp.c<8>(0x0f + 0x3fff);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 1, c8, 1, {uncomp, 0, 0});
    dp.update(i, 0x0f, 0x0f + 0x001f, c8, 14, {uncomp, 0, 1});
    dp.update(i, 0x0f, 0x0f + 0x3fff, c8, 23, {uncomp, 0, 2});
    dp.update_matrix(i, ofs_tab, len_tab, c0, 1,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b8_l ret(0x0e);
  if (input.size() > 0) ret.write<b1>(false);

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      if (li == 0) {
        ret.write<b1>(true);
      } else {
        ret.write<bnh, bnl>({4, 6}, {4, 0x0f});
        if (li == 1) {
          ret.write<b1, bnl>(true, {5, cmd.len - 0x0f});
        } else {
          ret.write<b1, bnl>(false, {14, cmd.len - 0x0f});
        }
      }
      ret.write<b8ln>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const auto wr = [&ret](const vrange& tp, size_t prefix_bits, size_t v) {
        const size_t x = tp.val + (v - tp.min);
        ret.write<bnh>({prefix_bits, x >> (tp.bitlen - prefix_bits)});
        ret.write<bnl>({tp.bitlen - prefix_bits, x});
      };
      ret.write<b1>(false);
      wr(len_tab[li], len_prefix[li], cmd.len);
      wr(ofs_tab[oi], ofs_prefix[oi], adr - cmd.lz_ofs());
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  std::reverse(ret.out.begin() + 0x0e, ret.out.end());
  ret.write<d16b>(0x0f);

  write32b(ret.out, 0, 0x43724d21); // CrM!
  write16b(ret.out, 4, 0); // unknown
  write32b(ret.out, 6, input.size());
  write32b(ret.out, 10, ret.size() - 14);
  assert(adr == input.size());
  assert(dp.optimal_cost() + (input.size() > 0 ? 1 : 0) + (14 + 2) * 8 == ret.bit_length());

  return ret.out;
}

} // namespace sfc_comp
