#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> rob_northen_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);
  enum Tag {
    uncomp, uncompl, lz
  };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag != lz) return true;
      return len_no == rhs.len_no && ofs_no == rhs.ofs_no; // ?
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  static constexpr size_t len_tab[][3] = {
    {0x0002, 0x0002,  2}, // 10
    {0x0003, 0x0006,  3}, // 110
    {0x0004, 0x0000,  3}, // 0_0
    {0x0006, 0x0002,  4}, // 0_1_
    {0x0009, 0x0007, 11}, // 111[]
    {0x0100, 0x0000,  0}
  };

  static constexpr size_t ofs_tab[][3] = {
    {0x0001, 0x0000,  9}, // 0[]
    {0x0101, 0x0006, 11}, // 110[]
    {0x0201, 0x0008, 12}, // 100_[]
    {0x0401, 0x0015, 13}, // 1_1_1[]
    {0x0801, 0x0028, 14}, // 1_1_0_[]
    {0x1001, 0x0000,  0}
  };

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_k<4>(i, 12, 72, Linear<8, 9>(), {uncompl, 0, 0});

    const auto cost = dp[i].cost;
    for (size_t k = 0; k < 5; ++k) {
      size_t oi = 4 - k;
      const size_t min_ofs = ofs_tab[oi][0];
      if (min_ofs > i) continue;
      const size_t max_ofs = ofs_tab[oi + 1][0] - 1;
      const size_t ofs_bitsize = ofs_tab[oi][2];
      auto res_lz = lz_helper.find_best(i, max_ofs);
      if (res_lz.len <= 1) break;
      if ((i - res_lz.ofs) < min_ofs) continue;
      for (size_t li = (oi == 0) ? 0 : 1; li < 5; ++li) {
        const size_t min_len = len_tab[li][0];
        if (min_len > res_lz.len) break;
        const size_t max_len = len_tab[li + 1][0] - 1;
        size_t len_bitsize = len_tab[li][2];
        size_t total_cost = cost + 1 + ofs_bitsize + len_bitsize;
        if (oi == 0 && li == 0) --total_cost;
        dp.update_lz(i, min_len, max_len, res_lz,
                     Constant<0>(), {lz, oi, li}, total_cost);
      }
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  for (size_t i = 0; i < 18; ++i) ret.write<d8>(0);
  ret.write<b1h, b1h>(false, false);

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1h, d8>(false, input[adr]);
    } break;
    case uncompl: {
      ret.write<b8hn_h, d8n>({9, 0x170 + (cmd.len - 12) / 4}, {cmd.len, &input[adr]});
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const size_t li = cmd.type.len_no;
      const size_t oi = cmd.type.ofs_no;

      ret.write<b1h>(true);
      const size_t len_bits = len_tab[li][2];
      if (li == 0) {
        ret.write<b8hn_h, d8>({len_bits, len_tab[li][1]},  d - 1);
      } else {
        const size_t min_len = len_tab[li][0];
        size_t v = cmd.len - min_len;
        if (li == 2) {
          v <<= 1;
        } else if (li == 3) {
          v += (v & ~1);
        }
        if (len_bits < 8) {
          ret.write<b8hn_h>({len_bits, v | len_tab[li][1]});
        } else {
          ret.write<b8hn_h>({len_bits - 8, len_tab[li][1]});
          ret.write<d8>(cmd.len - (min_len - 1));
        }
        size_t delta = d - ofs_tab[oi][0];
        size_t ov = delta >> 8;
        if (oi == 3) {
          ov <<= 1;
          ov += ov & ~3;
        } else if (oi == 4) {
          ov += ov & ~1;
          ov += ov & ~7;
        }
        ret.write<b8hn_h>({ofs_tab[oi][2] - 8, ov | ofs_tab[oi][1]});
        ret.write<d8>(delta & 0xff);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b8hn_h, d8, b1h>({4, 0x0f}, 0, false);
  assert(((18*8 + 2 + dp.total_cost() + 4 + 8 + 1) + 7) / 8 == ret.size());

  write32b(ret.out, 0, 0x524e4302);
  write32b(ret.out, 4, input.size());
  write32b(ret.out, 8, ret.size() - 18);
  write16b(ret.out, 12, utility::crc16(input, 0, input.size()));
  write16b(ret.out, 14, utility::crc16(ret.out, 18, ret.size() - 18));
  ret.out[0x10] = 0; // leeway
  ret.out[0x11] = 0; // chunks

  return ret.out;
}

} // namespace sfc_comp
