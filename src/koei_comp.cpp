#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> koei_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);
  enum Tag {
    uncomp, lz
  };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return true;
      return len_no == rhs.len_no && ofs_no == rhs.ofs_no; // ?
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  static constexpr size_t len_tab[][3] = {
    {0x0002, 0x0001,  1}, // 1
    {0x0003, 0x0002,  3}, // 01_
    {0x0005, 0x0004,  5}, // 001__
    {0x0009, 0x0008,  7}, // 0001___
    {0x0011, 0x0010,  9}, // 00001____
    {0x0021, 0x0020, 11}, // 000001_____
    {0x0041, 0x0040, 13}, // 0000001______
    {0x0081, 0x0000, 14}, // 0000000_______
    {0x0100, 0x0000,  0}
  };

  static constexpr size_t ofs_tab[][3] = {
    {0x0001, 0x0000,  6}, // 0000__
    {0x0005, 0x0008,  7}, // 00010__
    {0x0009, 0x0018,  8}, // 00011___
    {0x0021, 0x0060,  9}, // 0011_____
    {0x0081, 0x0180, 10}, // 011_______
    {0x0101, 0x0400, 11}, // 100________
    {0x0201, 0x0a00, 12}, // 101_________
    {0x0401, 0x1800, 13}, // 110__________
    {0x0801, 0x3800, 14}, // 111___________
    {0x1001, 0x0000, 0}
  };

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    const auto cost = dp[i].cost;
    for (size_t k = 0; k < 9; ++k) {
      size_t oi = 8 - k;
      const size_t min_ofs = ofs_tab[oi][0];
      if (min_ofs > i) continue;
      const size_t max_ofs = ofs_tab[oi + 1][0] - 1;
      const size_t ofs_bitsize = ofs_tab[oi][2];
      auto res_lz = lz_helper.find_best(i, max_ofs);
      if (res_lz.len <= 1) break;
      if ((i - res_lz.ofs) < min_ofs) continue;
      for (size_t li = 0; li < 8; ++li) {
        const size_t min_len = len_tab[li][0];
        if (min_len > res_lz.len) break;
        const size_t max_len = len_tab[li + 1][0] - 1;
        const size_t len_bitsize = len_tab[li][2];
        dp.update_lz(i, min_len, max_len, res_lz,
                     Constant<0>(), {lz, oi, li}, cost + 1 + ofs_bitsize + len_bitsize);
      }
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret; ret.write<d16>(0);

  size_t curr16 = 0, next16 = 0, bit16_pos = 0;
  auto write_b16 = [&](size_t input_bits, size_t v) {
    while (input_bits > 0) {
      --input_bits;
      if (!bit16_pos) {
        bit16_pos = 16;
        curr16 = next16; next16 = ret.size(); ret.write<d16>(0);
      }
      --bit16_pos;
      if ((v >> input_bits) & 1) {
        ret.out[curr16 + (bit16_pos >> 3)] |= 1 << (bit16_pos & 7);
      }
    }
  };

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1h, d8>(true, input[adr]); break;
    }
    case lz: {
      const size_t d = adr - cmd.lz_ofs, l = cmd.len;
      const size_t li = cmd.type.len_no;
      const size_t oi = cmd.type.ofs_no;
      ret.write<b1h>(false);
      write_b16(len_tab[li][2], len_tab[li][1] + (l - len_tab[li][0]));
      write_b16(ofs_tab[oi][2], ofs_tab[oi][1] + (d - ofs_tab[oi][0]));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1h>(false);
  write_b16(14, 0x007f);

  if (ret.size() == next16 + 2) {
    ret.out.resize(ret.size() - 2);
  }
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
