#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> rob_northen_comp_2_core(
    std::span<const uint8_t> input,
    const size_t header_size, std::span<const vrange> ofs_tab) {
  check_size(input.size(), 0, 0x800000);
  enum Tag { uncomp, uncompl, lz, lz2 };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag != lz) return true;
      return li == rhs.li;
    }
    Tag tag;
    size_t oi, li;
  };

  static constexpr std::array<vrange, 4> len_tab = {
                                        // 10 (len == 2)
    vrange(0x0003, 0x0003,  3, 0x0006), // 110
    vrange(0x0004, 0x0005,  3, 0x0000), // 0_0
    vrange(0x0006, 0x0008,  4, 0x0002), // 0_1_
    vrange(0x0009, 0x00ff, 11, 0x0007)  // 111[]
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_k<4>(i, 12, 72, Linear<8, 9>(), {uncompl, 0, 0});
    const auto res_lz2 = lz_helper.find_best(i, ofs_tab[0].max);
    dp.update_lz(i, 2, 2, res_lz2, Constant<11>(), {lz2, 0, 0});
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> CompType { return {lz, oi, li}; },
      1
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret(header_size);
  ret.write<b1, b1>(false, false);

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, d8>(false, input[adr]);
    } break;
    case uncompl: {
      ret.write<bnh, d8n>({9, 0x170 + (cmd.len - 12) / 4}, {cmd.len, &input[adr]});
    } break;
    case lz2: {
      ret.write<b1, bnh, d8>(true, {2, 2}, d - 1);
    } break;
    case lz: {
      const size_t li = cmd.type.li;
      const auto& l = len_tab[li];
      ret.write<b1>(true);
      size_t v = cmd.len - l.min;
      if (li == 1) {
        v <<= 1;
      } else if (li == 2) {
        v += (v & ~1);
      }
      if (l.bitlen < 8) {
        ret.write<bnh>({l.bitlen, v | l.val});
      } else {
        ret.write<bnh>({l.bitlen - 8, l.val});
        ret.write<d8>(cmd.len - (l.min - 1));
      }

      const size_t oi = cmd.type.oi;
      const auto& o = ofs_tab[oi];
      size_t delta = d - o.min;
      size_t ov = delta >> 8;
      if (oi == 3) {
        ov <<= 1;
        ov += ov & ~3;
      } else if (oi == 4) {
        ov += ov & ~1;
        ov += ov & ~7;
      }
      ret.write<bnh>({o.bitlen - 8, ov | o.val});
      ret.write<d8>(delta & 0xff);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh, d8, b1>({4, 0x0f}, 0, false);
  assert(header_size * 8 + 2 + dp.total_cost() + 4 + 8 + 1 == ret.bit_length());
  return ret.out;
}

static constexpr std::array<vrange, 5> rnc2_offsets = {
  vrange(0x0001, 0x0100,  9, 0x0000), // 0[]
  vrange(0x0101, 0x0200, 11, 0x0006), // 110[]
  vrange(0x0201, 0x0400, 12, 0x0008), // 100_[]
  vrange(0x0401, 0x0800, 13, 0x0015), // 1_1_1[]
  vrange(0x0801, 0x1000, 14, 0x0028), // 1_1_0_[]
};

} // namespace

std::vector<uint8_t> rob_northen_comp_2(std::span<const uint8_t> input) {
  auto ret = rob_northen_comp_2_core(input, 0x12, rnc2_offsets);
  write32b(ret, 0, 0x524e4302);
  write32b(ret, 4, input.size());
  write32b(ret, 8, ret.size() - 18);
  write16b(ret, 12, utility::crc16(input, 0, input.size()));
  write16b(ret, 14, utility::crc16(ret, 18, ret.size() - 18));
  ret[0x10] = 0; // leeway
  ret[0x11] = 0; // chunks
  return ret;
}

std::vector<uint8_t> spirou_comp(std::span<const uint8_t> input) {
  auto ret = rob_northen_comp_2_core(input, 0x02, rnc2_offsets);
  write16b(ret, 0, 0x524e);
  return ret;
}

std::vector<uint8_t> smurfs_comp(std::span<const uint8_t> input) {
  auto ret = rob_northen_comp_2_core(input, 0x02, std::span(rnc2_offsets.begin(), 3));
  write16b(ret, 0, 0x524e);
  return ret;
}

} // namespace sfc_comp
