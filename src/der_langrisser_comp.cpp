#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> der_langrisser_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lzs, lzl };
  using tag = tag_l<method>;

  static constexpr auto lens_s = to_vranges({
    {1, 1, 0b0},
    {2, 2, 0b10},
    {3, 3, 0b110},
  }, 3);

  static constexpr auto lens_l = to_vranges({
    {2, 4, 0b1110},
    {3, 5, 0b11110},
    {4, 6, 0b111110},
    {5, 11, 0b1111110'0000},
  }, 20);

  static constexpr auto codewords = create_array<encode::codeword, 256>([](size_t i) -> encode::codeword {
    if (i < 0x10) return {7, 0x40 | i};
    if (i == 0x10) return {5, 0x14};
    if (i == 0x30) return {5, 0x15};
    if (i == 0x80) return {5, 0x16};
    if (i == 0xff) return {5, 0x17};
    return {10, 0x300 | i};
  });

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(lens_l.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, codewords[input[i]].bitlen, {uncomp, 0});
    const auto res_lz4 = lz_helper.find(i, 0x10, lens_s.front().min);
    dp.update(i, lens_s, res_lz4, c0, 5, [&](size_t li) -> tag { return {lzs, li}; });
    const auto res_lz8 = lz_helper.find(i, 0x100, lens_l.front().min);
    dp.update(i, lens_l, res_lz8, c0, 9, [&](size_t li) -> tag { return {lzl, li}; });
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, li] = cmd.type;
    const auto d = adr - cmd.lz_ofs();
    switch (tag) {
    case uncomp: {
      ret.write<bnh>(codewords[input[adr]]);
    } break;
    case lzs: case lzl: {
      const auto& l = (tag == lzs) ? lens_s[li] : lens_l[li];
      ret.write<b1>(false);
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      if (tag == lzs) ret.write<bnh>({4, d});
      else ret.write<bnh>({8, d});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1, bnh>(false, {7, 0b1111111});
  assert(adr == input.size());
  assert(dp.optimal_cost() + 8 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
