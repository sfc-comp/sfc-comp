#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> syndicate_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum Tag {
    uncomp, lz1, lz2, lz3
  };
  struct CompType {
    Tag tag;
    size_t from;
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return true;
      return from == rhs.from;
    }
  };
  static constexpr size_t min_ofs_bits = 1, max_ofs_bits = 16;
  static constexpr size_t min_len_bits = 2, max_len_bits = 16;
  size_t len_masks[17] = {};
  for (size_t i = min_len_bits; i <= max_len_bits; ++i) {
    len_masks[i] = 0x0003 << (i - min_len_bits);
  }
  size_t min_lens[max_ofs_bits + max_len_bits + 1] = {};
  for (size_t i = 0; i < 8; ++i) min_lens[i] = 1;
  for (size_t i = 8; i < 17; ++i) min_lens[i] = 2;
  for (size_t i = 17; i <= max_len_bits + max_ofs_bits; ++i) min_lens[i] = 3;

  size_t max_offsets[max_ofs_bits + 1] = {};
  for (size_t i = min_ofs_bits; i <= max_ofs_bits; ++i) max_offsets[i] = (1 << i) - 1;

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, max_ofs_bits + 1>> lz_memo(input.size());

  size_t longest_lz_len = 0, longest_lz_dist = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    size_t o = max_ofs_bits;
    auto res_lz = lz_helper.find_best_closest(i, max_offsets[o], input.size());
    longest_lz_len = std::max(longest_lz_len, res_lz.len);
    if (res_lz.len > 0) longest_lz_dist = std::max(longest_lz_dist, i - res_lz.ofs);
    lz_memo[i][o] = res_lz;
    for (o -= 1; o > 0; --o) {
      const size_t d = i - res_lz.ofs;
      if (res_lz.len > 0 && d > max_offsets[o]) {
        res_lz = lz_helper.find_best_closest(i, max_offsets[o], input.size());
      }
      lz_memo[i][o] = res_lz;
    }
    lz_helper.add_element(i);
  }

  std::vector<uint8_t> best;

  for (size_t ofs_bits = min_ofs_bits; ofs_bits <= max_ofs_bits; ++ofs_bits) {
    size_t len_bits_limit = min_len_bits;
    for (size_t len_bits = min_len_bits; len_bits < max_len_bits; ++len_bits) {
      const size_t max_len = min_lens[len_bits + ofs_bits] + ((size_t(1) << len_bits) - 1);
      if (longest_lz_len > max_len) len_bits_limit = len_bits + 1;
    }
    std::vector<sssp_solver<CompType>> dp(len_bits_limit + 1);
    for (size_t len_bits = min_len_bits; len_bits <= len_bits_limit; ++len_bits) {
      dp[len_bits] = sssp_solver<CompType>(input.size());
    }
    for (size_t i = 0; i < input.size(); ++i) {
      auto res_lz = lz_memo[i][ofs_bits];
      for (size_t len_bits = min_len_bits; len_bits <= len_bits_limit; ++len_bits) {
        dp[len_bits].update(i, 1, 1, Constant<9>(), {uncomp, len_bits});
        const size_t curr_cost = dp[len_bits][i].cost;
        const size_t min_len = min_lens[len_bits + ofs_bits];
        const size_t next_cost = curr_cost + (1 + ofs_bits + len_bits);
        const size_t hi = len_masks[len_bits], lo = (hi & -hi);
        assert(min_len > 0);
        dp[std::max(min_len_bits, len_bits - 1)].update_lz(
          i, min_len, min_len + (lo - 1), res_lz, Constant<0>(), {lz1, len_bits}, next_cost);
        dp[len_bits].update_lz(
          i, min_len + lo, min_len + hi - 1, res_lz, Constant<0>(), {lz2, len_bits}, next_cost);
        dp[std::min(len_bits_limit, len_bits + 1)].update_lz(
          i, min_len + hi, std::min<size_t>(0x10000, min_len + (size_t(1) << len_bits) - 1), res_lz,
          Constant<0>(), {lz3, len_bits}, next_cost);
      }
    }

    size_t last_len_bits = 0, min_cost = dp[min_len_bits].infinite_cost;
    for (size_t len_bits = min_len_bits; len_bits <= len_bits_limit; ++len_bits) {
      if (dp[len_bits].total_cost() < min_cost) {
        min_cost = dp[len_bits].total_cost();
        last_len_bits = len_bits;
      }
    }
    size_t curr_len_bits = last_len_bits;
    const auto commands = [&input, &dp](size_t& curr) {
      using command_type = sssp_solver<CompType>::vertex_type;
      std::vector<command_type> ret;
      ptrdiff_t adr = input.size();
      while (adr > 0) {
        auto cmd = dp[curr][adr];
        assert(cmd.len > 0);
        adr -= cmd.len;
        curr = cmd.type.from;
        ret.emplace_back(cmd);
      }
      assert(adr == 0);
      std::reverse(ret.begin(), ret.end());
      return ret;
    }(curr_len_bits);

    using namespace data_type;
    writer_b ret; ret.write<b8ln_l>({32, 0});
    ret.out[0] = curr_len_bits; // initial len bits
    ret.out[1] = ofs_bits;
    ret.out[2] = min_lens[curr_len_bits + ofs_bits] - 1;
    ret.out[3] = len_bits_limit;

    size_t adr = 0;
    for (const auto& cmd : commands) {
      switch (cmd.type.tag) {
      case uncomp: {
        ret.write<b1l, b8ln_l>(true, {8, input[adr]});
      } break;
      case lz1: case lz2: case lz3: {
        const size_t d = adr - cmd.lz_ofs;
        const size_t min_len = min_lens[ofs_bits + curr_len_bits];
        assert(cmd.len >= min_len);
        const size_t mask = len_masks[curr_len_bits];
        const size_t l = cmd.len - min_len;
        ret.write<b1l, b8ln_l, b8ln_l>(false, {ofs_bits, d}, {curr_len_bits, l});
        if (l & mask) {
          if ((l & mask) == mask && curr_len_bits < len_bits_limit) ++curr_len_bits;
        } else if (curr_len_bits > min_len_bits) {
          --curr_len_bits;
        }
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    ret.write<b1l, b8ln_l>(0, {ofs_bits, 0});
    assert(adr == input.size());
    assert((min_cost + (1 + ofs_bits) + 7) / 8 + 4 == ret.size());
    if (best.empty() || ret.size() < best.size()) {
      best = std::move(ret.out);
    }
    if (max_offsets[ofs_bits] >= longest_lz_dist) break;
  }
  return best;
}

} // namespace sfc_comp
