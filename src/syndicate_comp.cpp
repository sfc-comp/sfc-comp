#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> syndicate_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz };
  using tag = tag_l<method>;

  static constexpr size_t ob_min = 1, ob_max = 16;
  static constexpr size_t lb_min = 2, lb_max = 16;
  static constexpr auto len_masks = create_array<size_t, lb_max + 1>([&](size_t i) {
    return (i < lb_min) ? 0 : 3 << (i - lb_min);
  });
  static constexpr auto min_lens = create_array<size_t, lb_max + ob_max + 1>([&](size_t i) {
    return (i < 9) ? 1
         : (i < 17) ? 2 : 3;
  });
  static constexpr size_t lz_min_len = 1;
  static constexpr auto max_offsets = create_array<size_t, ob_max + 1>([&](size_t i) {
    return (i < ob_min) ? 0 : (1 << i) - 1;
  });

  const auto [lz_memo, longest_lz_len, longest_lz_dist] = [&] {
    size_t longest_lz_len = 0, longest_lz_dist = 0;
    std::vector<std::array<encode::lz_data, max_offsets.size()>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, max_offsets, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find(i, max_ofs, lz_min_len);
      });
      if (const auto lz = ret[i].back(); lz.len >= lz_min_len) {
        longest_lz_len = std::max(longest_lz_len, lz.len);
        longest_lz_dist = std::max(longest_lz_dist, i - lz.ofs);
      }
      lz_helper.add_element(i);
    }
    return std::make_tuple(std::move(ret), longest_lz_len, longest_lz_dist);
  }();

  std::vector<uint8_t> best;

  for (size_t ob = ob_min; ob <= ob_max; ++ob) {
    size_t lb_lim = lb_min;
    for (size_t lb = lb_min; lb < lb_max; ++lb) {
      const size_t max_len = min_lens[lb + ob] + ((size_t(1) << lb) - 1);
      if (longest_lz_len > max_len) lb_lim = lb + 1;
    }
    std::vector<solver<tag>> dp(lb_lim + 1);
    for (size_t lb = lb_min; lb <= lb_lim; ++lb) dp[lb] = solver<tag>(input.size());
    std::vector<decltype(dp[0].c<0>(0))> c0s(lb_lim + 1);
    for (size_t lb = lb_min; lb <= lb_lim; ++lb) c0s[lb] = dp[lb].c<0>(((1 << lb) - 1) + 3);

    for (size_t i = input.size(); i-- > 0; ) {
      for (size_t lb = lb_min; lb <= lb_lim; ++lb) {
        dp[lb].update(i, 1, 9, {uncomp, lb});
        const size_t lbs = std::max(lb_min, lb - 1), lbl = std::min(lb_lim, lb + 1);
        const size_t mn = min_lens[lb + ob];
        const size_t hi = len_masks[lb], lo = (hi & -hi), mx = (hi | (lo - 1));
        const auto res_lz = lz_memo[i][ob];
        const size_t c = 1 + lb + ob;
        dp[lb].update(i, mn,      mn + lo - 1, res_lz, c0s[lbs], c, {lz, lbs});
        dp[lb].update(i, mn + lo, mn + hi - 1, res_lz, c0s[lb],  c, {lz, lb});
        dp[lb].update(i, mn + hi, mn + mx    , res_lz, c0s[lbl], c, {lz, lbl});
      }
      for (size_t lb = lb_min; lb <= lb_lim; ++lb) c0s[lb].update(i);
    }

    const auto [best_lb, min_cost] = [&] {
      size_t best_lb = -1, best = dp[lb_min].infinite_cost;
      for (size_t lb = lb_min; lb <= lb_lim; ++lb) {
        if (dp[lb].optimal_cost() < best) best = dp[lb].optimal_cost(), best_lb = lb;
      }
      return std::make_pair(best_lb, best);
    }();

    using namespace data_type;
    writer_b8_l ret(4);

    size_t adr = 0;
    for (size_t lb = best_lb; adr < input.size(); ) {
      const auto& cmd = dp[lb][adr];
      const auto [tag, li] = cmd.type;
      switch (tag) {
      case uncomp: {
        ret.write<b1, bnl>(true, {8, input[adr]});
      } break;
      case lz: {
        const size_t d = adr - cmd.lz_ofs();
        const size_t min_len = min_lens[ob + lb];
        assert(min_len <= cmd.len && cmd.len <= min_len + ((1 << lb) - 1));
        const size_t mask = len_masks[lb];
        const size_t l = cmd.len - min_len;
        ret.write<b1, bnl, bnl>(false, {ob, d}, {lb, l});
        if (l & mask) {
          if ((l & mask) == mask && lb < lb_lim) ++lb;
        } else if (lb > lb_min) --lb;
      } break;
      default: assert(0);
      }
      assert(lb == li);
      adr += cmd.len;
    }
    ret[0] = best_lb; // initial len bits
    ret[1] = ob;
    ret[2] = min_lens[best_lb + ob] - 1;
    ret[3] = lb_lim;
    ret.write<b1, bnl>(0, {ob, 0});
    assert(adr == input.size());
    assert(min_cost + (1 + ob) + 4 * 8 == ret.bit_length());
    if (best.empty() || ret.size() < best.size()) best = std::move(ret.out);
    if (max_offsets[ob] >= longest_lz_dist) break;
  }
  return best;
}

} // namespace sfc_comp
