#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> oscar_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  enum Tag { uncomp, lz };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return true;
      return li == rhs.li;
    }
    Tag tag;
    size_t oi, li;
  };

  static constexpr auto ofs_tab = create_array<vrange, 15>([](size_t k) {
    return vrange(1 << k, (2 << k) - 1, 4 + k, (k + 1) << k);
  });
  static constexpr auto ulen_tab = create_array<vrange, 8>([](size_t k) {
    return (k == 0) ? vrange(1, 1, 3, 0)
                    : vrange((1 << (k - 1)) + 1, (1 << k), 3 + (k - 1), k << (k - 1));
  });
  static constexpr size_t min_len = 3;
  static constexpr auto len_tab = create_array<vrange, ulen_tab.size()>([](size_t k) {
    return (k == 0) ? vrange(min_len, min_len, 3, 0)
                    : vrange((1 << (k - 1)) + min_len, (1 << k) + (min_len - 1), 3 + (k - 1), k << (k - 1));
  });

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    u_helper.update(i, dp[i].cost);
    for (size_t k = 0; k < ulen_tab.size(); ++k) {
      const auto u = u_helper.find(i + 1, ulen_tab[k].min, ulen_tab[k].max);
      dp.update_u(i + 1, u.len, {uncomp, 0, k}, u.cost + 4 + ulen_tab[k].bitlen);
    }
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> CompType { return {lz, oi, li}; },
      0
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret(4);
  writer raw;
  size_t adr = 0, counter = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    const size_t li = cmd.type.li;
    switch (cmd.type.tag) {
    case uncomp: {
      const auto& l = ulen_tab[li];
      ret.write<bnh, bnh>({4, 0}, {l.bitlen, l.val + (cmd.len - l.min)});
      raw.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[cmd.type.oi];
      ret.write<bnh>({o.bitlen, o.val + (d - o.min)});
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
    ++counter;
  }
  assert(adr == input.size());
  assert(dp.total_cost() + 32 == raw.size() * 8 + ret.bit_length());
  write16b(ret.out, 0, ret.size() - 4);
  write16b(ret.out, 2, counter);
  std::copy(raw.out.begin(), raw.out.end(), std::back_inserter(ret.out));
  return ret.out;
}

} // namespace sfc_comp
