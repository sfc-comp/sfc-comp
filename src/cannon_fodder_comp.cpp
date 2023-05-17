#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> cannon_fodder_comp(std::span<const uint8_t> input) {
  static constexpr size_t shift = 16;
  check_size(input.size(), 1, (1 << shift) + 1);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr size_t min_oi = 1, max_oi = 14;
  static constexpr size_t min_li = 4, max_li = 14;
  static constexpr auto len_masks = create_array<size_t, max_li + 1>([&](size_t i) {
    return (i < min_li) ? 0 : 0x000c << (i - min_li);
  });
  static constexpr size_t lz_min_len = 1;
  static constexpr auto len_vals =
      create_array<std::array<size_t, 4>, max_li + 1>([&](size_t li) -> std::array<size_t, 4> {
    if (li < min_li) return {};
    const size_t lo = (len_masks[li] & -len_masks[li]), hi = len_masks[li];
    return {lz_min_len, lo + lz_min_len, hi + lz_min_len, (1 << li) + lz_min_len};
  });
  static constexpr auto max_offsets = create_array<size_t, max_oi + 1>([&](size_t i) {
    return 1 << i;
  });

  const auto split = [](size_t oi, size_t adr) -> std::pair<size_t, size_t> {
    if (oi == max_oi) return {0, adr};
    uint64_t mask = low_bits_mask(oi);
    return {(adr >> oi) & 1, ((adr >> 1) & ~mask) | (adr & mask)};
  };
  const auto merge = [](size_t oi, size_t b, size_t adr) -> size_t {
    if (oi == max_oi) return adr;
    uint64_t mask = low_bits_mask(oi);
    return ((adr << 1) - (adr & mask)) | (b << oi);
  };
  const auto lowers = [&split](size_t oi, size_t i) -> std::pair<size_t, size_t> {
    if (oi == max_oi) return {i, 0};
    const auto mask = low_bits_mask(oi);
    auto [b, adr] = split(oi, i);
    if (b == 0) return {adr, adr & ~mask};
    return {adr + 1 + (~adr & mask), adr};
  };
  const auto uppers = [&split](size_t oi, size_t i) -> std::pair<size_t, size_t> {
    if (oi == max_oi) return {i, -1};
    const auto mask = low_bits_mask(oi);
    auto [b, adr] = split(oi, i);
    if (b == 0) return {adr, (adr & ~mask) - 1};
    return {adr | mask, adr};
  };

  const auto [lz_memo, longest_lz_len] = [&] {
    size_t longest_lz_len = lz_min_len;
    std::vector<std::array<encode::lz_data, max_offsets.size()>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, max_offsets, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find_closest(i, max_ofs, lz_min_len, input.size());
      });
      if (const auto lz = ret[i].back(); lz.len >= lz_min_len) {
        longest_lz_len = std::max(longest_lz_len, lz.len);
      }
      lz_helper.add_element(i);
    }
    return std::make_tuple(std::move(ret), longest_lz_len);
  }();

  const size_t li_limit = [&] {
    size_t ret = min_li;
    for (size_t li = 0; li < max_li; ++li) {
      if ((longest_lz_len - 1) >= len_masks[li]) ret = li + 1;
    }
    return ret;
  }();
  const size_t oi_limit = std::min(max_oi, std::bit_width(input.size()));

  using solver_type = solver<tag>;
  std::array<std::array<solver_type, max_oi + 1>, max_li + 1> dp;
  std::array<std::array<std::array<cost_window<0>, 2>, max_oi + 1>, max_li + 1> c0;

  for (size_t li = min_li; li <= li_limit; ++li) {
    for (size_t oi = min_oi; oi <= oi_limit; ++oi) {
      dp[li][oi] = solver_type(input.size());
      const auto [bit, size] = split(oi, input.size());
      c0[li][oi][bit] = cost_window<0>(size, len_vals[li].back() - 1);
      c0[li][oi][1 - bit] = cost_window<0>(input.size() - size, len_vals[li].back() - 1, -1);
    }
  }

  for (size_t i = input.size(); i-- > 0; ) {
    for (size_t li = min_li; li <= li_limit; ++li) {
      const auto nlis = std::to_array({std::max(min_li, li - 1), li, std::min(max_li, li + 1)});
      for (size_t oi = min_oi; oi <= oi_limit; ++oi) {
        const bool reachable = (oi == max_oi) || !((i >> oi) & 1) || (oi > min_oi && ((i >> (oi - 1)) & 1));
        const auto [b, adr] = split(oi, i);
        if (reachable) {
          const auto [nb, nadr] = split(oi, i + 1);
          dp[li][oi].update_c(i, 1, c0[li][oi][nb][nadr] + 9, {uncomp, std::min(oi + nb, max_oi), li});
          const auto res_lz = lz_memo[i][oi];
          if (res_lz.len >= lz_min_len) {
            const auto [len0, len1] = uppers(oi, i + res_lz.len);
            for (size_t k = 0; k < 3; ++k) {
              const auto [fr0, fr1] = lowers(oi, i + len_vals[li][k]);
              const auto [to0, to1] = uppers(oi, i + len_vals[li][k + 1] - 1);
              const size_t nli = nlis[k];
              if (~to0 && fr0 <= len0) {
                const auto e = c0[nli][oi][0].find(0, fr0, std::min(len0, to0));
                const auto l = merge(oi, 0, e.len) - i;
                dp[li][oi].update_c(i, l, e.cost + oi + li + 1, {lz, oi, nli}, res_lz.ofs);
              }
              if (~to1 && ~len1 && fr1 <= len1) {
                const size_t noi = std::min(oi + 1, max_oi); assert(noi <= oi_limit);
                const auto e = c0[nli][oi][1].find(0, fr1, std::min(len1, to1));
                const auto l = merge(oi, 1, e.len) - i;
                dp[li][oi].update_c(i, l, e.cost + oi + li + 1, {lz, noi, nli}, res_lz.ofs);
              }
            }
          }
          c0[li][oi][b].update(adr, dp[li][oi][i].cost);
          if (oi > min_oi && ((i >> (oi - 1)) & 1)) {
            const auto [lb, ladr] = split(oi - 1, i);
            c0[li][oi - 1][lb].update(ladr, dp[li][oi][i].cost);
          }
        } else {
          c0[li][oi][b].update(adr, solver<tag>::infinite_cost);
        }
      }
    }
  }

  using namespace data_type;
  writer_b8_h ret(4);
  size_t adr = 0;
  for (size_t li = min_li, oi = min_oi; adr < input.size(); ) {
    const auto& cmd = dp[li][oi][adr];
    const auto [tag, noi, nli] = cmd.type;
    assert(cmd.len > 0);
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, bnh>(false, {8, input[adr]});
    } break;
    case lz: {
      ret.write<b1, bnh, bnh>(true, {li, cmd.len - lz_min_len}, {oi, (adr - cmd.lz_ofs()) - 1});
      const size_t mask = len_masks[li];
      const size_t l = cmd.len - lz_min_len;
      if (l & mask) {
        if ((l & mask) == mask && li < max_li) ++li;
      } else if (li > min_li) --li;
    } break;
    default: assert(0);
    }
    adr += cmd.len;
    // This strange behavior complicates the algorithm.
    if (oi < max_oi && adr & (1 << oi)) oi += 1;
    assert(oi == noi && li == nli);
  }
  write32b(ret.out, 0, input.size());
  assert(adr == input.size());
  assert(dp[min_li][min_oi].optimal_cost() + 8 * 4 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
