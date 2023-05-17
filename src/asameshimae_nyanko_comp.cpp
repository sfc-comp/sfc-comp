#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "image.hpp"

namespace sfc_comp {

namespace {

struct nyanko_cost {
  constexpr nyanko_cost(size_t cost = 0, size_t lz_count = 0) : value(cost), lz_count(lz_count) {}
  constexpr auto operator <=> (const nyanko_cost& rhs) const = default;
  constexpr nyanko_cost operator + (size_t c) const { return nyanko_cost(value + c, lz_count); }
  constexpr nyanko_cost operator - (size_t c) const { return nyanko_cost(value - c, lz_count); }
  size_t value;
  size_t lz_count;
};

} // namespace

template <>
struct cost_traits<nyanko_cost> {
  static constexpr nyanko_cost infinity() { return nyanko_cost(cost_traits<size_t>::infinity()); }
  static constexpr nyanko_cost unspecified() { return nyanko_cost(cost_traits<size_t>::unspecified()); }
};

std::vector<uint8_t> asameshimae_nyanko_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);

  enum method { uncomp, lz, none, lz_cnt };
  using tag = tag_ol<method>;

  static constexpr size_t lz_min_len = 1;
  static constexpr auto ofs_tab = create_array<vrange, 16>([](size_t k) {
    return vrange(1 << k, (2 << k) - 1, 2 * k + 1, ((1 << k) - 1) << (k + 1));
  });
  static constexpr auto ulen_tab = create_array<vrange, 16>([](size_t k) {
    return vrange(1 << k, (2 << k) - 1, 2 * k + 1, ((1 << k) - 1) << (k + 1));
  });
  static constexpr auto len_tab = create_array<vrange, 16>([](size_t k) {
    return vrange((1 << k) + lz_min_len - 1, (2 << k) + lz_min_len - 2, 2 * k + 1, ((1 << k) - 1) << (k + 1));
  });
  static constexpr auto lz_lv = create_array<vrange, 16>([](size_t k) {
    return vrange(1 << k, (2 << k) - 1, 2 * k + 1, ((1 << k) - 1) << (k + 1));
  });
  static constexpr size_t lz_max_len = len_tab.back().max;

  if (input.size() > ulen_tab.back().max) {
    throw std::runtime_error("This algorithm may not be able to compress the given data.");
  }

  const auto lz_memo = [&] {
    std::vector<std::array<encode::lz_data, ofs_tab.size()>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, ofs_tab, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find_closest(i, max_ofs, lz_min_len, lz_max_len);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  const size_t lv_max = ilog2(input.size());

  solver<tag> dpu(input.size()), dpl(input.size());
  using solver_type = solver<tag, nyanko_cost>;
  std::vector<solver_type> dp(lv_max + 1, solver_type(input.size(), -1));

  auto c8 = dpl.c<8>(ulen_tab.back().max);
  using window_type = cost_window<0, 1, std::greater<size_t>, nyanko_cost>;
  std::vector<window_type> c0s(lv_max + 1, window_type(input.size(), len_tab.back().max, -1));

  size_t lv_lim = 0;
  for (size_t i = input.size(); ; ) {
    c0s[0].update(i, nyanko_cost(dpu[i].cost + lz_lv[0].bitlen, 1));
    if (i-- == 0) break;
    dpu.update(i, ulen_tab, c8, 0, [&](size_t li) -> tag { return {uncomp, 0, li}; });
    for (size_t lv = 0; lv <= lv_lim; ++lv) {
      dp[lv].update_matrix(i, ofs_tab, len_tab, c0s[lv], 0,
        [&](size_t oi) { return lz_memo[i][oi]; },
        [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
      );
    }
    for (size_t lv = 0, lv_e = lv_lim; lv <= lv_e; ++lv) {
      auto c = dp[lv][i].cost;
      if (c == dp[lv].infinite_cost) continue;
      dpl.update_c(i, 0, c.value, {lz_cnt, 0, lv}, c.lz_count);
      size_t nlv = lv;
      if (++c.lz_count > lz_lv[lv].max) {
        nlv += 1; assert(nlv <= lv_max);
        c.value += lz_lv[nlv].bitlen - lz_lv[lv].bitlen;
        if (nlv > lv_lim) lv_lim = nlv;
      }
      if (c < c0s[nlv][i]) c0s[nlv].update(i, c);
    }
    c8.update(i);
  }

  using namespace data_type;
  writer ret(5); writer_b8_h flags;

  size_t adr = 0;
  for (size_t curr = 0; adr < input.size(); ) {
    if (curr < 2) {
      const auto& cmd = (curr == 0) ? dpu[adr] : dpl[adr];
      const auto [tag, oi, li] = cmd.type;
      switch (tag) {
      case uncomp: {
        const auto& l = ulen_tab[li];
        flags.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
        ret.write<d8n>({cmd.len, &input[adr]});
        curr = 1;
      } break;
      case lz_cnt: {
        const auto& lv = lz_lv[li];
        assert(lv.min <= cmd.arg && cmd.arg <= lv.max);
        flags.write<bnh>({lv.bitlen, lv.val + (cmd.arg - lv.min)});
        curr = li + 2;
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    } else {
      const auto& cmd = dp[curr - 2][adr];
      const auto [tag, oi, li] = cmd.type;
      switch (tag) {
      case lz: {
        const size_t d = adr - cmd.lz_ofs();
        const auto& l = len_tab[li];
        const auto& o = ofs_tab[oi];
        flags.write<bnh>({o.bitlen, o.val + (d - o.min)});
        flags.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
        if (cmd.cost.lz_count == lz_lv[curr - 2].min) {
          if (--curr == 1) curr = 0;
        }
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
  }
  assert(adr == input.size());
  assert(dpu.optimal_cost() + 5 * 8 == ret.size() * 8 + flags.bit_length());
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 5);
  ret[4] = 0;
  ret.extend(flags);

  return ret.out;
}

std::vector<uint8_t> asameshimae_nyanko_4bpp_comp(std::span<const uint8_t> input) {
  check_divisibility(input.size(), 0x20);
  auto ret = asameshimae_nyanko_comp(snes4bpp::to_indexed16_h_8_1(input));
  ret[4] = 1;
  return ret;
}

} // namespace sfc_comp
