#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> cannon_fodder_comp(std::span<const uint8_t> input) {
  static constexpr size_t shift = 16;
  check_size(input.size(), 1, (1 << shift) + 1);

  enum method { none, uncomp, lzs, lzm, lzl };
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
      for (auto& lz : ret[i]) lz.ofs = lz.ofs << shift | i;
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

  using solver_type = sssp_solver<tag>;
  std::array<std::array<std::array<solver_type, 2>, max_oi + 1>, max_li + 1> dp;
  for (size_t li = min_li; li <= li_limit; ++li) {
    for (size_t oi = min_oi; oi <= oi_limit; ++oi) {
      const auto [bit, size] = split(oi, input.size());
      dp[li][oi][bit] = solver_type(size, -1);
      dp[li][oi][1 - bit] = solver_type(input.size() - size, -1);
    }
  }
  dp[min_li][min_oi][0][0].cost = 0;

  // [TODO] Too slow for *simple* inputs.
  const auto infinite_cost = solver_type::infinite_cost;
  const auto tags = std::to_array<method>({lzs, lzm, lzl});
  for (size_t i = 0; i < input.size(); ++i) {
    for (size_t li = min_li; li <= li_limit; ++li) {
      const auto nlis = std::to_array({std::max(min_li, li - 1), li, std::min(max_li, li + 1)});
      for (size_t oi = min_oi; oi <= oi_limit; ++oi) {
        const auto [b, adr] = split(oi, i);
        const auto cost = dp[li][oi][b][adr].cost;
        if (cost == infinite_cost) continue;
        size_t coi = oi + (b > 0); assert(coi <= max_oi);
        const auto [nb, nadr] = split(coi, i + 1);
        dp[li][coi][nb].update(0, nadr, nadr, Constant<9>(), {uncomp, oi, li}, cost, i);

        const auto lz = lz_memo[i][coi];
        if (lz.len < lz_min_len) continue;
        const auto [len0, len1] = uppers(coi, i + lz.len);
        const encode::lz_data lz0 = {lz.ofs, len0};
        const encode::lz_data lz1 = {lz.ofs, len1};
        for (size_t k = 0; k < 3; ++k) {
          const auto [fr0, fr1] = lowers(coi, i + len_vals[li][k]);
          const auto [to0, to1] = uppers(coi, i + len_vals[li][k + 1] - 1);
          if (~to0 && fr0 <= lz0.len) {
            dp[nlis[k]][coi][0].update_lz(0, fr0, to0, lz0, Constant<0>(), {tags[k], oi, li}, cost + coi + li + 1);
          }
          if (~to1 && ~lz1.len && fr1 <= lz1.len) {
            dp[nlis[k]][coi][1].update_lz(0, fr1, to1, lz1, Constant<0>(), {tags[k], oi, li}, cost + coi + li + 1);
          }
        }
      }
    }
  }

  const auto [commands, min_cost] = [&] {
    using command_type = solver_type::vertex_type;
    size_t best_cost = infinite_cost;
    size_t best_li = 0, best_oi = 0;
    for (size_t li = min_li; li <= li_limit; ++li) {
      for (size_t oi = min_oi; oi <= oi_limit; ++oi) {
        const auto [b, adr] = split(oi, input.size());
        if (dp[li][oi][b][adr].cost < best_cost) {
          best_cost = dp[li][oi][b][adr].cost;
          best_li = li, best_oi = oi;
        }
      }
    }
    std::vector<command_type> ret;
    size_t li = best_li, oi = best_oi;
    auto [b, adr] = split(oi, input.size());
    while (!(b == 0 && adr == 0)) {
      auto cmd = dp[li][oi][b][adr];
      size_t poi = cmd.type.oi;
      size_t nadr = cmd.val() & low_bits_mask(shift);
      cmd.lz_ofs = cmd.val() >> shift;
      cmd.len = merge(oi, b, adr) - nadr;
      std::tie(b, adr) = split(poi, nadr);
      li = cmd.type.li; oi = poi;
      ret.emplace_back(cmd);
    }
    assert(li == min_li && oi == min_oi && b == 0 && adr == 0);
    std::ranges::reverse(ret);
    return std::make_pair(std::move(ret), best_cost);
  }();

  using namespace data_type;
  writer_b8_h ret(4);
  size_t adr = 0;
  size_t curr_li = min_li, curr_oi = min_oi;
  for (const auto& cmd : commands) {
    switch (cmd.type.tag) {
    case none: break;
    case uncomp: {
      ret.write<b1, bnh>(false, {8, input[adr]});
    } break;
    case lzs: case lzm: case lzl: {
      ret.write<b1, bnh, bnh>(true, {curr_li, cmd.len - 1}, {curr_oi, (adr - cmd.lz_ofs) - 1});
      const size_t mask = len_masks[curr_li];
      const size_t l = cmd.len - 1;
      if (l & mask) {
        if ((l & mask) == mask && curr_li < max_li) ++curr_li;
      } else if (curr_li > min_li) --curr_li;
    } break;
    default: assert(0);
    }
    adr += cmd.len;
    // This strange behavior complicates the algorithm.
    if (curr_oi < max_oi && adr & (1 << curr_oi)) curr_oi += 1;
  }
  write32b(ret.out, 0, input.size());
  assert(adr == input.size());
  assert(min_cost + 8 * 4 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
