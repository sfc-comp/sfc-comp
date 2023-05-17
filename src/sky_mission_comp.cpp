#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> sky_mission_comp_core(std::span<const uint8_t> input, const bool mixed) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  7, 0b00'00000},
    {0x0021,  9, 0b01'0000000},
    {0x00a1, 11, 0b10'000000000},
    {0x02a1, 12, 0b11'0000000000}
  }, 0x06a0);

  static constexpr auto len_tab = to_vranges({
    {0x0002,  1, 0b0},                    // 0
    {0x0003,  3, 0b1'00 + 1},             // 1__
    {0x0006,  7, 0b100'0000 + 1},         // 100____
    {0x0015, 15, 0b1000000'00000000 + 1}  // 1000000________
  }, 0x0113);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.c<0>(len_tab.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, {uncomp, 0, 0});
    dp.update_matrix(i, ofs_tab, len_tab, c0, 1,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  writer raw;

  if (!mixed) ret.write<d16>(0);

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      ret.write<b1>(false);
      if (mixed) {
        ret.write<bnh>({8, input[adr]});
      } else {
        raw.write<d8>(input[adr]);
      }
    } break;
    case lz: {
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[oi];
      ret.write<b1>(true);
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      ret.write<bnh>({o.bitlen, o.val + ((adr - cmd.lz_ofs()) - o.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1>(true);
  ret.write<bnh>({15, 1 << 14});
  assert(adr == input.size());
  assert(dp.optimal_cost() + 1 + 15 + (mixed ? 0 : 2) * 8 == 8 * raw.size() + ret.bit_length());

  if (!mixed) {
    if (ret.size() >= 0x10000) throw std::runtime_error("This algorithm cannot compress the given data.");
    write16(ret.out, 0, ret.size());
    ret.extend(raw);
  }
  return ret.out;
}

} // namespace

std::vector<uint8_t> sky_mission_comp(std::span<const uint8_t> input) {
  return sky_mission_comp_core(input, true);
}

std::vector<uint8_t> riddick_bowe_boxing_comp(std::span<const uint8_t> input) {
  return sky_mission_comp_core(input, false);
}

} // namespace sfc_comp
