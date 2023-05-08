#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> soul_and_sword_comp(std::span<const uint8_t> input) {
  static constexpr size_t skipped_size = 2;

  check_size(input.size(), skipped_size, 0x10000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto uncomp_len_tab = to_vranges({
    {0x0001,  3, 0b000},
    {0x0004,  4, 0b0110},
    {0x0005,  6, 0b011101},
    {0x0006,  8, 0b0111101'0},
    {0x0008, 10, 0b01111101'00},
    {0x000c, 12, 0b011111101'000},
    {0x0014, 14, 0b0111111101'0000},
    {0x0024, 16, 0b01111111101'00000},
    {0x0044, 18, 0b011111111101'000000},
    {0x0084, 20, 0b0111111111101'0000000},
  }, 0x0100);

  static constexpr auto lz_len_tab = to_vranges({
    {0x0002,  3, 0b000},
    {0x0009,  4, 0b1110},
    {0x000a,  6, 0b111101},
    {0x000b,  8, 0b1111101'0},
    {0x000d, 10, 0b11111101'00},
    {0x0011, 12, 0b111111101'000},
    {0x0019, 14, 0b1111111101'0000},
    {0x0029, 16, 0b11111111101'00000},
    {0x0049, 18, 0b111111111101'000000},
    {0x0089, 20, 0b1111111111101'0000000},
  }, 0x00ff);

  static constexpr auto ofs_tab = to_vranges({
    {0x0001, 0, 0} // dynamically changed
  }, 0x0fff);

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size(), skipped_size);

  for (size_t i = 0; i < skipped_size; ++i) lz_helper.add_element(i);

  for (size_t i = skipped_size; i < input.size(); ++i) {
    const auto cost = dp[i].cost;

    for (size_t k = 0; k < uncomp_len_tab.size(); ++k) {
      const size_t min_len = uncomp_len_tab[k].min;
      if (i + min_len > input.size()) break;
      const size_t max_len = uncomp_len_tab[k].max;
      const size_t len_bitsize = uncomp_len_tab[k].bitlen;
      const size_t total_cost = cost + len_bitsize;
      dp.update(i, min_len, max_len, Linear<8, 0>(), {uncomp, 0, k}, total_cost);
    }
    const size_t ofs_bitsize = std::min<size_t>(12, 1 + ilog2(i));
    dp.update_lz_matrix(i, ofs_tab, lz_len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t, size_t li) -> tag { return {lz, ofs_bitsize, li}; },
      1 + ofs_bitsize
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret(4);
  writer raws; raws.write<d8n>({skipped_size, &input[0]});

  size_t adr = skipped_size;
  for (const auto& cmd : dp.commands(skipped_size)) {
    switch (cmd.type.tag) {
    case uncomp: {
      const auto& l = uncomp_len_tab[cmd.type.li];
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      raws.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const auto& l = lz_len_tab[cmd.type.li];
      ret.write<b1>(true);
      assert(d < (size_t(1) << cmd.type.oi));
      ret.write<bnh>({cmd.type.oi, d});
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 4);
  assert(dp.total_cost() + 8 * (4 + skipped_size) == raws.size() * 8 + ret.bit_length());
  ret.extend(raws);

  return ret.out;
}

} // namespace sfc_comp
