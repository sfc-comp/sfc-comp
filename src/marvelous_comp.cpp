#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> marvelous_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum method { uncomp = 0, rle = 1, rle16 = 2, inc = 3, lz = 4 };
  using tag = tag_l<method>;
  static constexpr auto lens = to_vranges({{0x0001, 1, 0}, {0x0021, 2, 0}, {0x0401, 3, 0}}, 0x10000);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x10000);
  auto c1 = dp.c<1>(0x10000);

  size_t rlen = 0, rlen16 = 0, rleni = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, lens, c1, 0, [&](size_t li) -> tag { return {uncomp, li}; });
    rlen = encode::run_length_r(input, i, rlen);
    dp.update(i, lens, rlen, c0, 1, [&](size_t li) -> tag { return {rle, li}; });
    rlen16 = encode::run_length16_r(input, i, rlen16);
    dp.update(i, lens, rlen16, c0, 2, [&](size_t li) -> tag { return {rle16, li}; });
    rleni = encode::run_length_r(input, i, rleni, 1);
    dp.update(i, lens, rleni, c0, 1, [&](size_t li) -> tag { return {inc, li}; });
    dp.update(i, lens, lz_helper.find(i, 0x10000, 3), c0, 2, [&](size_t li) -> tag { return {lz, li}; });

    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, li] = cmd.type;
    if (li == 0) ret.write<d8>(tag << 5 | (cmd.len - 1));
    else if (li == 1) ret.write<d16b>(0xe000 | (tag << 10) | (cmd.len - 1));
    else ret.write<d24b>(0xc00000 | (tag << 18) | (cmd.len - 1));
    switch (tag) {
    case uncomp: ret.write<d8n>({cmd.len, &input[adr]}); break;
    case rle: ret.write<d8>(input[adr]); break;
    case rle16: ret.write<d8, d8>(input[adr], input[adr + 1]); break;
    case inc: ret.write<d8>(input[adr]); break;
    case lz: ret.write<d16b>(cmd.lz_ofs()); break;
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
