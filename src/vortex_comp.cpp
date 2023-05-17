#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> vortex_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0x800000);

  enum method { none, uncomp, lz2, lz3, lz };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = to_vranges({
             // 00 (len == 2)
             // 01 (len == 3)
    {  4,  2, 0b10},
    {  5,  4, 0b110'0},
    {  7,  6, 0b1110'00},
    { 11, 12, 0b1111'00000000 + 0x0b}
  }, 255);

  static constexpr auto ofs_tab = to_vranges({
    {0x0001, 10, 0b11'00000000 + 0x0001},
    {0x0100, 14, 0b10'000000000000 + 0x0100},
    {0x1000, 17, 0b0'0000000000000000 + 0x1000}
  }, 0xffff);

  static constexpr auto ulen_tab = to_vranges({
    {0x0001,  3, 0b000 + 1},
    {0x0007,  8, 0b1110'0000},
    {0x0017, 14, 0b1111'0000000000 + 0x17}
  }, 0x03ff);

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input, true);
  solver<tag> dp0(input.size(), -1); auto c0_0 = dp0.c<0>(len_tab.back().max);
  solver<tag> dp1(input.size()); auto c8_1 = dp1.c<8>(ulen_tab.back().max);

  for (size_t i = input.size(); ; ) {
    dp0.update_c(i, 0, dp1[i].cost + 3, {none, 0, 0});
    c0_0.update(i);

    if (i-- == 0) break;
    lz_helper.reset(i);

    dp0.update(i, ulen_tab, c8_1, 0, [&](size_t li) -> tag { return {uncomp, 0, li}; });

    const auto res_lzs = lz_helper.find(i, 0xff, 2);
    dp1.update(i, 2, 2, res_lzs, c0_0, 10, {lz2, 0, 0});
    dp1.update(i, 3, 3, res_lzs, c0_0, 11, {lz3, 0, 0});
    dp1.update(i, 3, 3, lz_helper.find(i, 0x3fff, 3),  c0_0, 17, {lz3, 1, 0});
    dp1.update_matrix(i, ofs_tab, len_tab, c0_0, 0,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c8_1.update(i);
  }

  if (dp0.optimal_cost() == dp0.infinite_cost) {
    // For example, this happens when the input is a de Bruijn sequence.
    throw std::runtime_error("This algorithm cannot compress the given data.");
  }

  using namespace data_type;
  writer_b8_l ret(4); ret.write<b1>(false);
  size_t adr = 0;
  for (size_t curr = 0; adr <= input.size(); curr ^= 1) {
    if (adr == input.size() && curr == 1) break;
    const auto& cmd = (curr == 0) ? dp0[adr] : dp1[adr];
    const auto [tag, oi, li] = cmd.type;
    const size_t d = adr - cmd.lz_ofs();
    switch (tag) {
      case none: {
        ret.write<bnh>({3, 0b000});
      } break;

      case uncomp: {
        const auto& u = ulen_tab[li];
        ret.write<bnh>({u.bitlen, u.val + (cmd.len - u.min)});
        ret.write<b8hn>({cmd.len, &input[adr]});
      } break;

      case lz2: {
        ret.write<bnh, bnh>({2, 0b00}, {8, d});
      } break;

      case lz3: {
        ret.write<bnh>({2, 0b01});
        if (d < 0x100) ret.write<bnh, bnh>({1, 0b1}, {8, d});
        else ret.write<bnh, bnh>({1, 0b0}, {14, d});
      } break;

      case lz: {
        const auto& l = len_tab[li];
        const auto& o = ofs_tab[oi];
        ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
        ret.write<bnh>({o.bitlen, o.val + (d - o.min)});
      } break;

      default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp0.optimal_cost() + 33 == ret.bit_length());
  write32(ret.out, 0, input.size());

  const size_t s = std::min<size_t>(8, ret.size());
  uint64_t v = 1;
  for (size_t i = s - 1; i >= 4; --i) v = v << 8 | ret.out[i];
  v >>= 1;
  for (size_t i = 4; i < s; ++i) ret.out[i] = v & 0xff, v >>= 8;

  std::ranges::reverse(ret.out);
  return ret.out;
}

} // namespace sfc_comp
