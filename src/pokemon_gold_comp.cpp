#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> pokemon_gold_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x8000);

  enum method { uncomp = 0, rle = 1, rle16 = 2, rle0 = 3, lz = 4, lzh = 5, lzv = 6 };
  using tag = tag_ol<method>;

  static constexpr auto lens = to_vranges({{0x0001, 1, 0}, {0x0021, 2, 0}}, 0x0400);
  static constexpr auto offsets = to_vranges({{0x0001, 1, 0}, {0x0081, 2, 0}}, 0x8000);

  lz_helper_kirby lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x400);
  auto c1 = dp.c<1>(0x400);

  size_t rlen = 0, rlen16 = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, lens, c1, 0, [&](size_t li) -> tag { return {uncomp, 0, li}; });
    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, lens, rlen, c0, 0, [&](size_t li) -> tag { return {rle0, 0, li}; });
    } else {
      dp.update(i, lens, rlen, c0, 1, [&](size_t li) -> tag { return {rle, 0, li}; });
    }
    rlen16 = encode::run_length16_r(input, i, rlen16);
    dp.update(i, lens, rlen16, c0, 2, [&](size_t li) -> tag { return {rle16, 0, li}; });

    dp.update_matrix(i, offsets, lens, c0, 0,
      [&](size_t oi) { return lz_helper.find(i, offsets[oi].max, 2); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    dp.update_matrix(i, offsets, lens, c0, 0,
      [&](size_t oi) { return lz_helper.find_h(i, offsets[oi].max, 2); },
      [&](size_t oi, size_t li) -> tag { return {lzh, oi, li}; }
    );
    dp.update_matrix(i, offsets, lens, c0, 0,
      [&](size_t oi) { return lz_helper.find_v(i, offsets[oi].max, 2); },
      [&](size_t oi, size_t li) -> tag { return {lzv, oi, li}; }
    );

    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    if (li == 0) ret.write<d8>(tag << 5 | (cmd.len - 1));
    else ret.write<d16b>(0xe000 | (tag << 10) | (cmd.len - 1));
    switch (tag) {
    case uncomp: ret.write<d8n>({cmd.len, &input[adr]}); break;
    case rle: ret.write<d8>(input[adr]); break;
    case rle16: ret.write<d8, d8>(input[adr], input[adr + 1]); break;
    case rle0: break;
    case lz: case lzh: case lzv: {
      if (oi == 0) ret.write<d8>(0x80 | ((adr - cmd.lz_ofs()) - 1));
      else ret.write<d16b>(cmd.lz_ofs());
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.optimal_cost() == ret.size());
  assert(adr == input.size());
  ret.write<d8>(0xff);
  return ret.out;
}

} // namespace sfc_comp
