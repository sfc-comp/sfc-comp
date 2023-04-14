#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> diet_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  enum Tag {
    uncomp, lz2, lz
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
    {0x0003, 0x0001,  1}, // 1
    {0x0004, 0x0001,  2}, // 01
    {0x0005, 0x0001,  3}, // 001
    {0x0006, 0x0001,  4}, // 0001
    {0x0007, 0x0002,  6}, // 00001_
    {0x0009, 0x0000,  9}, // 000000___
    {0x0011, 0x0001, 14}, // 000001<1B>
    {0x0111, 0x0000,  0}  // * end *
  };

  static constexpr size_t ofs_tab[][3] = {
    {0x0001, 0x0001, 10}, // _1
    {0x0201, 0x0001, 11}, // _01
    {0x0401, 0x0001, 13}, // _00_1
    {0x0801, 0x0001, 15}, // _00_0_1
    {0x1001, 0x0000, 16}, // _00_0_0_
    {0x2001, 0x0000,  0}  // * end *
  };

  encode::lz_data res_lz2 = {}, res_lz2s = {};
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});

    if (res_lz2.len >= 2) {
      if (res_lz2s.len >= 2) {
        dp.update_lz(i, 2, 2, res_lz2s, Constant<11>(), {lz2, 0, 0});
      } else {
        dp.update_lz(i, 2, 2, res_lz2, Constant<14>(), {lz2, 1, 0});
      }
    }

    const auto cost = dp[i].cost;
    for (size_t k = 0; k < 5; ++k) {
      size_t oi = 4 - k;
      const size_t min_ofs = ofs_tab[oi][0];
      if (min_ofs > i) continue;
      const size_t max_ofs = ofs_tab[oi + 1][0] - 1;
      const size_t ofs_bitsize = ofs_tab[oi][2];
      auto res_lz = lz_helper.find_best(i, max_ofs);
      if (res_lz.len <= 2) break;
      if ((i - res_lz.ofs) < min_ofs) continue;
      for (size_t li = 0; li < 7; ++li) {
        const size_t min_len = len_tab[li][0];
        if (min_len > res_lz.len) break;
        const size_t max_len = len_tab[li + 1][0] - 1;
        size_t len_bitsize = len_tab[li][2];
        size_t total_cost = cost + 2 + ofs_bitsize + len_bitsize;
        dp.update_lz(i, min_len, max_len, res_lz,
                     Constant<0>(), {lz, oi, li}, total_cost);
      }
    }

    if (i + 1 < input.size()) {
      res_lz2 = lz_helper.find_best(i + 1, 0x900);
      if (res_lz2.len >= 2) {
        res_lz2s = lz_helper.find_best(i + 1, 0x100);
      }
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  for (size_t i = 0; i < 17; ++i) ret.write<d8>(0);
  ret.write<d16>(0);

  size_t curr16 = 0x11, bit16_pos = 0;
  auto write_b16 = [&](size_t input_bits, size_t v) {
    while (input_bits > 0) {
      --input_bits;
      if ((v >> input_bits) & 1) {
        ret.out[curr16 + (bit16_pos >> 3)] |= 1 << (bit16_pos & 7);
      }
      ++bit16_pos;
      if (bit16_pos == 16) {
        bit16_pos = 0;
        curr16 = ret.size();
        ret.write<d16>(0);
      }
    }
  };

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      write_b16(1, 1);
      ret.write<d8>(input[adr]);
    } break;
    case lz2: {
      const size_t d = adr - cmd.lz_ofs;
      write_b16(2, 0);
      ret.write<d8>((0x900 - d) & 0xff);
      if (cmd.type.ofs_no == 0) {
        write_b16(1, 0);
      } else {
        write_b16(4, 0x08 | (0x900 - d) >> 8);
      }
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const size_t li = cmd.type.len_no;
      const size_t oi = cmd.type.ofs_no;

      write_b16(2, 1);
      ret.write<d8>(-d & 0xff);
      size_t v = (ofs_tab[oi + 1][0] - 1 - d) >> 8;
      if (oi == 0) {
        v <<= 1;
      } else if (oi == 1) {
        v <<= 2;
      } else if (oi == 2) {
        v <<= 1;
        v += (v & 0x04) * 3;
      } else if (oi == 3) {
        v <<= 1;
        v += v & ~3;
        v += (v & 0x10) * 3;
      } else {
        v += v & ~1;
        v += v & ~7;
        v += (v & 0x20) * 3;
      }
      write_b16(ofs_tab[oi][2] - 8, v | ofs_tab[oi][1]);

      const size_t len_bits = len_tab[li][2];
      const size_t len_v = cmd.len - len_tab[li][0];
      if (len_bits < 10) {
        write_b16(len_tab[li][2], len_v | len_tab[li][1]);
      } else {
        write_b16(len_tab[li][2] - 8, len_tab[li][1]);
        ret.write<d8>(len_v);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write_b16(2, 0);
  ret.write<d8>(0xff);
  write_b16(1, 0);
  if (ret.size() == curr16 + 2 && bit16_pos == 0) {
    ret.out.resize(ret.size() - 2);
  }
  assert(adr == input.size());

  std::array<uint8_t, 9> header = {
    0xb4, 0x4c, 0xcd, 0x21, 0x9d, 0x89, 0x64, 0x6c, 0x7a
  };
  for (size_t i = 0; i < header.size(); ++i) ret.out[i] = header[i];
  ret[9] = ((ret.size() - 0x11) >> 16) & 0x0f;
  write16(ret.out, 10, ret.size() - 0x11);
  write16(ret.out, 12, utility::crc16(ret.out, 0x11, ret.size() - 0x11));
  ret.out[14] = ((input.size() >> 16) & 0x3f) << 2;
  write16(ret.out, 15, input.size());

  return ret.out;
}

} // namespace sfc_comp
