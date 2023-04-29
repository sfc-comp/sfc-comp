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
      return li == rhs.li;
    }
    Tag tag;
    size_t oi, li;
  };

  const std::array<vrange, 9> ofs_tab = {
    vrange(0x0001, 0x0004,  6, 0x0000), // 0000__
    vrange(0x0005, 0x0008,  7, 0x0008), // 00010__
    vrange(0x0009, 0x0020,  8, 0x0018), // 00011___
    vrange(0x0021, 0x0080,  9, 0x0060), // 0011_____
    vrange(0x0081, 0x0100, 10, 0x0180), // 011_______
    vrange(0x0101, 0x0200, 11, 0x0400), // 100________
    vrange(0x0201, 0x0400, 12, 0x0a00), // 101_________
    vrange(0x0401, 0x0800, 13, 0x1800), // 110__________
    vrange(0x0801, 0x1000, 14, 0x3800)  // 111___________
  };

  const std::array<vrange, 8> len_tab = {
    vrange(0x0002, 0x0002,  1, 0x0001), // 1
    vrange(0x0003, 0x0004,  3, 0x0002), // 01_
    vrange(0x0005, 0x0008,  5, 0x0004), // 001__
    vrange(0x0009, 0x0010,  7, 0x0008), // 0001___
    vrange(0x0011, 0x0020,  9, 0x0010), // 00001____
    vrange(0x0021, 0x0040, 11, 0x0020), // 000001_____
    vrange(0x0041, 0x0080, 13, 0x0040), // 0000001______
    vrange(0x0081, 0x00ff, 14, 0x0000)  // 0000000_______
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> CompType { return {lz, oi, li}; },
      1
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret(2);

  size_t curr16 = 0, next16 = 0, bit = 0;
  auto write_b16 = [&](size_t bitlen, size_t v) {
    while (bitlen > 0) {
      --bitlen;
      if (!bit) {
        bit = 16;
        curr16 = next16; next16 = ret.size(); ret.write<d16>(0);
      }
      --bit;
      if ((v >> bitlen) & 1) {
        ret.out[curr16 + (bit >> 3)] |= 1 << (bit & 7);
      }
    }
  };

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, d8>(true, input[adr]); break;
    }
    case lz: {
      const auto& l = len_tab[cmd.type.li];
      const auto& o = ofs_tab[cmd.type.oi];
      ret.write<b1>(false);
      write_b16(l.bitlen, l.val + (cmd.len - l.min));
      write_b16(o.bitlen, o.val + ((adr - cmd.lz_ofs) - o.min));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1>(false);
  write_b16(14, 0x007f);
  if (ret.size() == next16 + 2) {
    ret.out.resize(ret.size() - 2);
  }
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
