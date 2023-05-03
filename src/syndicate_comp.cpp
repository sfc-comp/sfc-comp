#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> syndicate_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz1, lz2, lz3 };
  using tag = tag_l<method>;

  static constexpr size_t ofs_min_bits = 1, ofs_max_bits = 16;
  static constexpr size_t len_min_bits = 2, len_max_bits = 16;
  static constexpr auto len_masks = create_array<size_t, len_max_bits + 1>([&](size_t i) {
    return (i < len_min_bits) ? 0 : 3 << (i - len_min_bits);
  });
  static constexpr auto min_lens = create_array<size_t, len_max_bits + ofs_max_bits + 1>([&](size_t i) {
    return (i < 9) ? 1
         : (i < 17) ? 2 : 3;
  });
  static constexpr size_t lz_min_len = 1;
  static constexpr auto max_offsets = create_array<size_t, ofs_max_bits + 1>([&](size_t i) {
    return (i < ofs_min_bits) ? 0 : (1 << i) - 1;
  });

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, ofs_max_bits + 1>> lz_memo(input.size());

  size_t longest_lz_len = 0, longest_lz_dist = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    for (ptrdiff_t oi = ofs_max_bits; oi >= 0; ) {
      auto res_lz = lz_helper.find_best_closest(i, max_offsets[oi], input.size());
      longest_lz_len = std::max(longest_lz_len, res_lz.len);
      if (res_lz.len >= lz_min_len) longest_lz_dist = std::max(longest_lz_dist, i - res_lz.ofs);
      else res_lz = {0, 0};
      do {
        lz_memo[i][oi--] = res_lz;
      } while (oi >= 0 && (res_lz.len < lz_min_len || (i - res_lz.ofs) <= max_offsets[oi]));
    }
    lz_helper.add_element(i);
  }

  std::vector<uint8_t> best;

  for (size_t ofs_bits = ofs_min_bits; ofs_bits <= ofs_max_bits; ++ofs_bits) {
    size_t len_bits_limit = len_min_bits;
    for (size_t len_bits = len_min_bits; len_bits < len_max_bits; ++len_bits) {
      const size_t max_len = min_lens[len_bits + ofs_bits] + ((size_t(1) << len_bits) - 1);
      if (longest_lz_len > max_len) len_bits_limit = len_bits + 1;
    }
    std::vector<sssp_solver<tag>> dp(len_bits_limit + 1);
    for (size_t len_bits = len_min_bits; len_bits <= len_bits_limit; ++len_bits) {
      dp[len_bits] = sssp_solver<tag>(input.size());
    }
    for (size_t i = 0; i < input.size(); ++i) {
      auto res_lz = lz_memo[i][ofs_bits];
      for (size_t len_bits = len_min_bits; len_bits <= len_bits_limit; ++len_bits) {
        dp[len_bits].update(i, 1, 1, Constant<9>(), {uncomp, len_bits});
        const size_t curr_cost = dp[len_bits][i].cost;
        const size_t min_len = min_lens[len_bits + ofs_bits];
        const size_t next_cost = curr_cost + (1 + ofs_bits + len_bits);
        const size_t hi = len_masks[len_bits], lo = (hi & -hi), mx = hi | (lo - 1);
        assert(min_len > 0);
        dp[std::max(len_min_bits, len_bits - 1)].update_lz(
          i, min_len, min_len + (lo - 1), res_lz, Constant<0>(), {lz1, len_bits}, next_cost);
        dp[len_bits].update_lz(
          i, min_len + lo, min_len + hi - 1, res_lz, Constant<0>(), {lz2, len_bits}, next_cost);
        dp[std::min(len_bits_limit, len_bits + 1)].update_lz(
          i, min_len + hi, std::min<size_t>(0x10000, min_len + mx), res_lz,
          Constant<0>(), {lz3, len_bits}, next_cost);
      }
    }

    size_t last_len_bits = 0, min_cost = dp[len_min_bits].infinite_cost;
    for (size_t len_bits = len_min_bits; len_bits <= len_bits_limit; ++len_bits) {
      if (dp[len_bits].total_cost() < min_cost) {
        min_cost = dp[len_bits].total_cost();
        last_len_bits = len_bits;
      }
    }
    size_t curr_len_bits = last_len_bits;
    const auto commands = [&input, &dp](size_t& curr) {
      using command_type = sssp_solver<tag>::vertex_type;
      std::vector<command_type> ret;
      ptrdiff_t adr = input.size();
      while (adr > 0) {
        auto cmd = dp[curr][adr];
        assert(cmd.len > 0);
        adr -= cmd.len;
        curr = cmd.type.li;
        ret.emplace_back(cmd);
      }
      assert(adr == 0);
      std::reverse(ret.begin(), ret.end());
      return ret;
    }(curr_len_bits);

    using namespace data_type;
    writer_b8_l ret(4);
    ret[0] = curr_len_bits; // initial len bits
    ret[1] = ofs_bits;
    ret[2] = min_lens[curr_len_bits + ofs_bits] - 1;
    ret[3] = len_bits_limit;

    size_t adr = 0;
    for (const auto& cmd : commands) {
      switch (cmd.type.tag) {
      case uncomp: {
        ret.write<b1, bnl>(true, {8, input[adr]});
      } break;
      case lz1: case lz2: case lz3: {
        const size_t d = adr - cmd.lz_ofs;
        const size_t min_len = min_lens[ofs_bits + curr_len_bits];
        assert(cmd.len >= min_len);
        const size_t mask = len_masks[curr_len_bits];
        const size_t l = cmd.len - min_len;
        ret.write<b1, bnl, bnl>(false, {ofs_bits, d}, {curr_len_bits, l});
        if (l & mask) {
          if ((l & mask) == mask && curr_len_bits < len_bits_limit) ++curr_len_bits;
        } else if (curr_len_bits > len_min_bits) {
          --curr_len_bits;
        }
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    ret.write<b1, bnl>(0, {ofs_bits, 0});
    assert(adr == input.size());
    assert(min_cost + (1 + ofs_bits) + 4 * 8 == ret.bit_length());
    if (best.empty() || ret.size() < best.size()) {
      best = std::move(ret.out);
    }
    if (max_offsets[ofs_bits] >= longest_lz_dist) break;
  }
  return best;
}

} // namespace sfc_comp
