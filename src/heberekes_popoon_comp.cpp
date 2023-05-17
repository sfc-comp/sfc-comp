#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> heberekes_popoon_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum tag { uncomp, uncompl, lz, lzpl, lzpo, lzpp };

  struct dist_len {
    size_t dist;
    size_t len;
  };

  static constexpr size_t max_dist = 0x100;
  static constexpr size_t max_len = 0x100;

  static constexpr size_t pre_dist_size = 0x0e;
  static constexpr size_t pre_len_size = 0x0f;

  lz_helper lz_helper(input);
  std::vector<encode::lz_data> lz_memo(input.size());

  size_t longest_lz_len = 0, longest_lz_dist = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    const auto res_lz = lz_helper.find(i, max_dist, 1);
    lz_memo[i] = res_lz;
    if (res_lz.len > 0) {
      longest_lz_dist = std::max(longest_lz_dist, i - res_lz.ofs);
      longest_lz_len = std::max(longest_lz_len, res_lz.len);
    }
    lz_helper.add_element(i);
  }
  longest_lz_len = std::min(longest_lz_len, max_len);
  longest_lz_dist = std::min(longest_lz_dist, max_dist);

  static constexpr auto coeffs = std::to_array<size_t>({256, 64, 32, 16, 8, 4, 2, 0});

  std::vector<dist_len> pre_sizes;
  for (size_t i = 0; i < coeffs.size(); ++i) {
    const auto f = [](size_t m, size_t v, size_t c) { return (v <= m) ? v : (m + (v - m) * c / 256); };
    const size_t len = f(pre_len_size, longest_lz_len, coeffs[i]);
    const size_t dist = f(pre_dist_size, longest_lz_dist, coeffs[i]);
    if (i > 0) {
      const auto last = pre_sizes.back();
      if (last.dist == dist && last.len == len) continue;
    }
    pre_sizes.push_back({dist, len});
  }

  std::vector<size_t> lens(longest_lz_len); std::iota(lens.begin(), lens.end(), 1);
  std::vector<size_t> dists(longest_lz_dist); std::iota(dists.begin(), dists.end(), 1);

  const size_t iter_total = pre_sizes.size();
  for (size_t iter = 0; iter < iter_total; ++iter) {
    solver<tag> dp(input.size());
    auto c0 = dp.c<0>(max_len);
    auto c1 = dp.c<1>(0x100);

    std::vector<size_t> res_lz(dists.size());
    for (size_t i = input.size(); i-- > 0; ) {
      dp.update(i, 0x10, 0x100, c1, 2, uncompl);
      dp.update(i, 0x01, 0x0f, c1, 1, uncomp);

      dp.update(i, 3, max_len, lz_memo[i], c0, 3, lz);
      dp.update(i, lens, lz_memo[i], constant<2>(), lzpl);

      if (dists.size() > 0) {
        for (size_t k = 0; k < dists.size(); ++k) {
          res_lz[k] = encode::lz_dist_r(input, i, dists[k], res_lz[k]);
        }
        const size_t best_k = std::max_element(res_lz.begin(), res_lz.end()) - res_lz.begin();
        const size_t best_len = res_lz[best_k];
        dp.update(i, lens, {best_k, best_len}, constant<1>(), lzpp);
        dp.update(i, 2, max_len, {best_k, best_len}, c0, 2, lzpo);
      }

      c0.update(i); c1.update(i);
    }

    const auto len_map = inverse_map<max_len + 1>(lens);
    if (iter + 1 < iter_total) {
      std::vector<size_t> dist_count(max_dist + 1);
      std::vector<size_t> len_count(max_len + 1);
      for (const auto& cmd : dp.optimal_path()) {
        if (cmd.type == lzpl) {
          len_count[cmd.len] += 1;
        } else if (cmd.type == lzpo) {
          dist_count[dists[cmd.arg]] += 1;
        } else if (cmd.type == lzpp) {
          len_count[cmd.len] += 1;
          dist_count[dists[cmd.arg]] += 1;
        }
      }
      const auto update = [](std::vector<size_t>& vals, size_t nsize, std::span<const size_t> counts) {
        std::partial_sort(vals.begin(), vals.begin() + nsize, vals.end(),
                          [&](size_t a, size_t b) { return counts[a] > counts[b]; });
        vals.resize(nsize);
        std::sort(vals.begin(), vals.end());
      };
      update(lens, pre_sizes[iter + 1].len, len_count);
      update(dists, pre_sizes[iter + 1].dist, dist_count);
    } else {
      using namespace data_type;
      writer ret(0x21);

      size_t adr = 0;
      for (const auto& cmd : dp.optimal_path()) {
        const size_t d = adr - cmd.lz_ofs();
        switch (cmd.type) {
        case lzpp: ret.write<d8>((cmd.arg << 4) | len_map[cmd.len]); break;
        case lzpo: ret.write<d8, d8>((cmd.arg << 4) | 0x0f, cmd.len & 0xff); break;
        case uncomp: ret.write<d8, d8n>(0xe0 | (cmd.len - 1), {cmd.len, &input[adr]}); break;
        case uncompl: ret.write<d8, d8, d8n>(0xef, cmd.len & 0xff, {cmd.len, &input[adr]}); break;
        case lzpl: ret.write<d8, d8>(0xf0 | len_map[cmd.len], d - 1); break;
        case lz: ret.write<d8, d8, d8>(0xff, d - 1, cmd.len & 0xff); break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      write32(ret.out, 0, input.size());
      for (size_t i = 0; i < dists.size(); ++i) ret[0x04 + i] = dists[i] - 1;
      for (size_t i = 0; i < lens.size(); ++i) ret[0x12 + i] = lens[i] & 0xff;
      assert(adr == input.size());
      assert(dp.optimal_cost() + 0x21 == ret.size());
      return ret.out;
    }
  }
  throw std::logic_error("iter_total == 0");
}

} // namespace sfc_comp
