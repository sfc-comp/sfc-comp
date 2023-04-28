#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sky_mission_comp_core(std::span<const uint8_t> input, const bool mixed) {
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

  static constexpr std::array<vrange, 4> ofs_tab = {
    vrange(0x0001, 0x0020,  7, 0x0000), // 00_____
    vrange(0x0021, 0x00a0,  9, 0x0080), // 01_______
    vrange(0x00a1, 0x02a0, 11, 0x0400), // 10_________
    vrange(0x02a1, 0x06a0, 12, 0x0c00)  // 11__________
  };

  static constexpr std::array<vrange, 4> len_tab = {
    vrange(0x0002, 0x0002,   1, 0x0000),     // 0
    vrange(0x0003, 0x0005,   3, 0x0004 + 1), // 1__
    vrange(0x0006, 0x0014,   7, 0x0040 + 1), // 100____
    vrange(0x0015, 0x0113,  15, 0x4000 + 1)  // 1000000________
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
  writer_b ret;
  writer raw;

  if (!mixed) ret.write<d16>(0);

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1h>(false);
      if (mixed) {
        ret.write<b8hn_h>({8, input[adr]});
      } else {
        raw.write<d8>(input[adr]);
      }
    } break;
    case lz: {
      const auto& l = len_tab[cmd.type.li];
      const auto& o = ofs_tab[cmd.type.oi];
      ret.write<b1h>(true);
      ret.write<b8hn_h>({l.bitlen, l.val + (cmd.len - l.min)});
      ret.write<b8hn_h>({o.bitlen, o.val + ((adr - cmd.lz_ofs) - o.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1h>(true);
  ret.write<b8hn_h>({15, 1 << 14});
  assert(adr == input.size());
  assert(dp.total_cost() + 1 + 15 + (mixed ? 0 : 2) * 8 == 8 * raw.size() + ret.bit_length());

  if (!mixed) {
    if (ret.size() >= 0x10000) throw std::runtime_error("This algorithm cannot compress the given data.");
    write16(ret.out, 0, ret.size());
    std::copy(raw.out.begin(), raw.out.end(), std::back_inserter(ret.out));
  }
  return ret.out;
}

std::vector<uint8_t> sky_mission_comp(std::span<const uint8_t> input) {
  return sky_mission_comp_core(input, true);
}

std::vector<uint8_t> riddick_bowe_boxing_comp(std::span<const uint8_t> input) {
  return sky_mission_comp_core(input, false);
}

} // namespace sfc_comp
