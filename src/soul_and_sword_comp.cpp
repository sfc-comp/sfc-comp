#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> soul_and_sword_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 2, 0x10000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto uncomp_len_tab = std::to_array<vrange>({
    vrange(0x0001, 0x0003,  3, 0b000),
    vrange(0x0004, 0x0004,  4, 0b0110),
    vrange(0x0005, 0x0005,  6, 0b011101),
    vrange(0x0006, 0x0007,  8, 0b0111101'0),
    vrange(0x0008, 0x000b, 10, 0b01111101'00),
    vrange(0x000c, 0x0013, 12, 0b011111101'000),
    vrange(0x0014, 0x0023, 14, 0b0111111101'0000),
    vrange(0x0024, 0x0043, 16, 0b01111111101'00000),
    vrange(0x0044, 0x0083, 18, 0b011111111101'000000),
    vrange(0x0084, 0x0100, 20, 0b0111111111101'0000000),
  });

  static constexpr auto lz_len_tab = std::to_array<vrange>({
    vrange(0x0002, 0x0008,  3, 0b000),
    vrange(0x0009, 0x0009,  4, 0b1110),
    vrange(0x000a, 0x000a,  6, 0b111101),
    vrange(0x000b, 0x000c,  8, 0b1111101'0),
    vrange(0x000d, 0x0010, 10, 0b11111101'00),
    vrange(0x0011, 0x0018, 12, 0b111111101'000),
    vrange(0x0019, 0x0028, 14, 0b1111111101'0000),
    vrange(0x0029, 0x0048, 16, 0b11111111101'00000),
    vrange(0x0049, 0x0088, 18, 0b111111111101'000000),
    vrange(0x0089, 0x00ff, 20, 0b1111111111101'0000000),
  });

  static constexpr auto ofs_tab = std::to_array<vrange>({
    vrange(0x0001, 0x0fff, 0, 0) // dynamically changed
  });

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size(), 2);

  for (size_t i = 0; i < 2; ++i) lz_helper.add_element(i);

  for (size_t i = 2; i < input.size(); ++i) {
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
  writer raws; raws.write<d8, d8>(input[0], input[1]);

  size_t adr = 2;
  for (const auto& cmd : dp.commands(2)) {
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
  assert(dp.total_cost() + 8 * (4 + 2) == raws.size() * 8 + ret.bit_length());
  std::copy(raws.out.begin(), raws.out.end(), std::back_inserter(ret.out));

  return ret.out;
}

} // namespace sfc_comp
