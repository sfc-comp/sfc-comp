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
      return len_no == rhs.len_no && ofs_no == rhs.ofs_no;
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  struct tup {
    size_t min_len;
    size_t val;
    size_t bits;
  };

  static constexpr std::array<tup, 5> len_tab = {
    tup({0x0002, 0x0000,      1}),  // 0
    tup({0x0003, 0x0004 + 1,  3}),  // 1__
    tup({0x0006, 0x0040 + 1,  7}),  // 100____
    tup({0x0015, 0x4000 + 1,  15}), // 1000000________
    tup({0x0114, 0x0000,      0})
  };

  static constexpr std::array<tup, 5> ofs_tab = {
    tup({0x0001, 0x0000, 7}),  // 00_____
    tup({0x0021, 0x0080, 9}),  // 01_______
    tup({0x00a1, 0x0400, 11}), // 10_________
    tup({0x02a1, 0x0c00, 12}), // 11__________
    tup({0x06a1, 0x0000, 0})
  };

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    const auto cost = dp[i].cost;
    for (ptrdiff_t oi = ofs_tab.size() - 2; oi >= 0; --oi) {
      const size_t min_ofs = ofs_tab[oi].min_len;
      if (min_ofs > i) continue;
      const size_t max_ofs = ofs_tab[oi + 1].min_len - 1;
      const size_t ofs_bitsize = ofs_tab[oi].bits;
      auto res_lz = lz_helper.find_best(i, max_ofs);
      if (res_lz.len <= 1) break;
      if ((i - res_lz.ofs) < min_ofs) continue;
      for (size_t li = 0; li < len_tab.size() - 1; ++li) {
        const size_t min_len = len_tab[li].min_len;
        if (min_len > res_lz.len) break;
        const size_t max_len = len_tab[li + 1].min_len - 1;
        const size_t len_bitsize = len_tab[li].bits;
        dp.update_lz(i, min_len, max_len, res_lz,
                     Constant<0>(), {lz, size_t(oi), li}, cost + 1 + ofs_bitsize + len_bitsize);
      }
    }
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
      const size_t d = adr - cmd.lz_ofs, l = cmd.len;
      const size_t li = cmd.type.len_no;
      const size_t oi = cmd.type.ofs_no;
      ret.write<b1h>(true);
      ret.write<b8hn_h>({len_tab[li].bits, len_tab[li].val + (l - len_tab[li].min_len)});
      ret.write<b8hn_h>({ofs_tab[oi].bits, ofs_tab[oi].val + (d - ofs_tab[oi].min_len)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1h>(true);
  ret.write<b8hn_h>({15, 1 << 14});
  assert(adr == input.size());
  assert((dp.total_cost() + 7 + 16) / 8 + (mixed ? 0 : 2) == raw.size() + ret.size());

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
