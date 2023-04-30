#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> diet_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  enum Tag { uncomp, lz2, lz };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag != lz) return true;
      return li == rhs.li;
    }
    Tag tag;
    size_t oi, li;
  };

  static constexpr auto ofs_tab = std::to_array<vrange>({
    vrange(0x0001, 0x0200, 10, 0b01),       // _1
    vrange(0x0201, 0x0400, 11, 0b001),      // _01
    vrange(0x0401, 0x0800, 13, 0b00001),    // _00_1
    vrange(0x0801, 0x1000, 15, 0b0000001),  // _00_0_1
    vrange(0x1001, 0x2000, 16, 0b00000000), // _00_0_0_
  });

  static constexpr auto len_tab = std::to_array<vrange>({
    vrange(0x0003, 0x0003,  1, 0b1),          // 1
    vrange(0x0004, 0x0004,  2, 0b01),         // 01
    vrange(0x0005, 0x0005,  3, 0b001),        // 001
    vrange(0x0006, 0x0006,  4, 0b0001),       // 0001
    vrange(0x0007, 0x0008,  6, 0b00001'0),    // 00001_
    vrange(0x0009, 0x0010,  9, 0b000000'000), // 000000___
    vrange(0x0011, 0x0110, 14, 0b000001)      // 000001<1B>
  });

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

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
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> CompType { return {lz, oi, li}; },
      2
    );
    // dist should be >= 2.
    if (i + 1 < input.size()) {
      res_lz2 = lz_helper.find_best(i + 1, 0x900);
      if (res_lz2.len >= 2) {
        res_lz2s = lz_helper.find_best(i + 1, 0x100);
      }
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16_hasty_l ret(17);
  ret.write<none>(none());

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, d8>(true, input[adr]);
    } break;
    case lz2: {
      const size_t d = adr - cmd.lz_ofs;
      ret.write<bnh, d8>({2, 0}, (0x900 - d) & 0xff);
      if (cmd.type.oi == 0) {
        ret.write<bnh>({1, 0});
      } else {
        ret.write<bnh>({4, 0x08 | (0x900 - d) >> 8});
      }
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const size_t oi = cmd.type.oi;
      const auto& o = ofs_tab[oi];

      ret.write<bnh, d8>({2, 1}, -d & 0xff);
      size_t v = (o.max - d) >> 8;
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
      ret.write<bnh>({o.bitlen - 8, v | o.val});

      const auto& l = len_tab[cmd.type.li];
      const size_t len_v = cmd.len - l.min;
      if (l.bitlen < 10) {
        ret.write<bnh>({l.bitlen, len_v | l.val});
      } else {
        ret.write<bnh, d8>({l.bitlen - 8, l.val}, len_v);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh, d8, bnh>({2, 0}, 0xff, {1, 0});
  ret.trim();

  assert(adr == input.size());
  assert(dp.total_cost() + 11 + 0x11 * 8 == ret.bit_length());

  static constexpr auto header = std::to_array<uint8_t>({
    0xb4, 0x4c, 0xcd, 0x21, 0x9d, 0x89, 0x64, 0x6c, 0x7a
  });
  std::copy(header.begin(), header.end(), ret.out.begin());
  ret[9] = ((ret.size() - 0x11) >> 16) & 0x0f;
  write16(ret.out, 10, ret.size() - 0x11);
  write16(ret.out, 12, utility::crc16(ret.out, 0x11, ret.size() - 0x11));
  ret.out[14] = ((input.size() >> 16) & 0x3f) << 2;
  write16(ret.out, 15, input.size());

  return ret.out;
}

} // namespace sfc_comp
