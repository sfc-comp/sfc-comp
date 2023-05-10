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

  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag, nyanko_cost> dpu(input.size()), dpl(input.size(), -1);
  std::vector<decltype(dpu)> dp(lv_max + 1, decltype(dpu)(input.size(), -1));

  for (size_t i = 0; i <= input.size(); ++i) {
    for (size_t lv = 0; lv <= lv_max; ++lv) {
      const auto cost = nyanko_cost(dp[lv][i].cost.value);
      if (cost == dpu.infinite_cost) continue;
      dpu.update_lz(i, 0, 0, {lv, 0}, Constant<0>(), {none, 0, 0}, cost);
    }
    if (i == input.size()) break;
    u_helper.update(i, dpu[i].cost.value);
    for (size_t k = 0; k < ulen_tab.size(); ++k) {
      const auto u = u_helper.find(i + 1, ulen_tab[k].min, ulen_tab[k].max);
      dpl.update_u(i + 1, u.len, {uncomp, 0, k}, u.cost + ulen_tab[k].bitlen);
    }
    const auto upd = [&](size_t nlv, size_t c, nyanko_cost cost) {
      dp[nlv].update_lz_matrix(i, ofs_tab, len_tab,
        [&](size_t oi) { return lz_memo[i][oi]; },
        [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
        c, cost
      );
    };
    const auto costl = dpl[i].cost;
    if (costl != dpu.infinite_cost) {
      upd(0, 1, nyanko_cost(costl.value, 1));
    }

    auto cutoff_v = dpu.infinite_cost.value;
    for (size_t lv = 0; lv <= lv_max; ++lv, cutoff_v += 2) {
      const auto cost = dp[lv][i].cost;
      if (cost == dpu.infinite_cost) continue;
      // This integer value is always even.
      // Therefore, this greedy algorithm produces an optimal output.
      const auto v = cost.value;
      if (v >= cutoff_v) continue;
      cutoff_v = v;
      const size_t nlv = (cost.lz_count == lz_lv[lv].max) ? lv + 1 : lv;
      if (nlv > lv_max) throw std::logic_error("lv_max too small.");
      upd(nlv, lz_lv[nlv].bitlen - lz_lv[lv].bitlen, nyanko_cost(cost.value, cost.lz_count + 1));
    }
  }

  const auto best_cost = std::min(dpu.optimal_cost(), dpl.optimal_cost()).value;
  const auto commands = [&] {
    using command_type = decltype(dpu)::vertex_type;
    std::vector<command_type> ret;
    size_t curr = (dpu.optimal_cost().value == best_cost) ? 0 : 1;
    ptrdiff_t adr = input.size();
    size_t lz_count = 0, lv = -1;
    while (adr > 0 || (adr == 0 && curr > 0)) {
      command_type cmd;
      if (curr == 0) {
        cmd = dpu[adr]; lv = cmd.val();
        lz_count = dp[lv][adr].cost.lz_count;
        curr = lv + 2;
      } else {
        if (curr == 1) {
          cmd = dpl[adr]; curr = 0;
          if (lz_count > 0) ret.emplace_back(-1, 0, lz_count, tag(lz_cnt, 0, lv));
        } else {
          cmd = dp[curr - 2][adr];
          if (cmd.cost.lz_count == lz_lv[curr - 2].min) curr -= 1;
        }
        assert(cmd.len > 0);
        adr -= cmd.len;
        ret.emplace_back(cmd);
      }
    }
    assert(adr == 0 && curr == 0);
    std::reverse(ret.begin(), ret.end());
    return ret;
  }();

  using namespace data_type;
  writer ret(5); writer_b8_h flags;

  size_t adr = 0;
  for (const auto& cmd : commands) {
    const size_t li = cmd.type.li;
    switch (cmd.type.tag) {
    case uncomp: {
      const auto& l = ulen_tab[li];
      flags.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz_cnt: {
      const auto& lv = lz_lv[li];
      assert(lv.min <= cmd.val() && cmd.val() <= lv.max);
      flags.write<bnh>({lv.bitlen, lv.val + (cmd.val() - lv.min)});
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[cmd.type.oi];
      flags.write<bnh>({o.bitlen, o.val + (d - o.min)});
      flags.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(best_cost + 5 * 8 == ret.size() * 8 + flags.bit_length());
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
