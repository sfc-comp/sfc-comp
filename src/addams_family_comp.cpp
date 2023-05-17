#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> addams_family_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz };
  using tag = tag_l<method>;

  static constexpr auto lens = to_vranges({
    {0x002,  1, 0b0},
    {0x003,  4, 0b1'000 + 1},
    {0x00a, 12, 0b1000'00000000 + 1}
  }, 0x108);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(9 + 0xff);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, {uncomp, 0});
    const auto res_lz = lz_helper.find(i, 0x400, 2);
    dp.update(i, lens, res_lz, c0, 1 + 10, [&](size_t li) -> tag { return {lz, li}; });
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      ret.write<b1, bnh>(false, {8, input[adr]});
    } break;
    case lz: {
      const auto& l = lens[li];
      ret.write<b1>(true);
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      ret.write<bnh>({10, adr - cmd.lz_ofs() - 1});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({13, 0x1800});
  assert(adr == input.size());
  assert(dp.optimal_cost() + 13 == ret.bit_length());
  return ret.out;
}

std::vector<uint8_t> jurassic_park_comp(std::span<const uint8_t> input) {
  auto ret = addams_family_comp(input);
  if (ret.size() & 1) ret.push_back(0);
  for (size_t i = 0; i < ret.size(); i += 2) std::swap(ret[i], ret[i + 1]);
  return ret;
}

} // namespace sfc_comp
