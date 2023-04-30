#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_donkey_kong_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  enum CompType {
    rle, uncomp, lz, lzl, pre16
  };

  static constexpr auto num_candidates = std::to_array<size_t>({
    2048, 512, 256, 128, 96, 80, 72, 68, 64
  });
  static constexpr size_t phase_total = num_candidates.size();

  auto candidate = utility::k_most_freq_u16(input, num_candidates[0]);
  std::vector<int64_t> pre(0x10000, -1);

  lz_helper lz_helper(input);
  std::vector<encode::lz_data> lz_memo(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    lz_memo[i] = lz_helper.find_best(i, 0xffff);
    lz_helper.add_element(i);
  }

  for (size_t phase = 0; phase < phase_total; ++phase) {
    for (size_t i = 0; i < candidate.size(); ++i) pre[candidate[i]] = i;

    sssp_solver<CompType> dp(input.size());

    size_t rlen = 0;
    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 0x3f, Linear<1, 1>(), uncomp);
      rlen = encode::run_length(input, i, rlen);
      dp.update(i, 1, 0x3f, rlen, Constant<2>(), rle);
      const auto res_lz = lz_memo[i];
      dp.update_lz(i, 1, 0x3f, res_lz, Constant<3>(), lz);
      dp.update_lz(i, 0x100, 0x100, res_lz, Constant<3>(), lzl);
      if (i + 1 < input.size()) {
        int16_t ind = pre[read16(input, i)];
        if (ind >= 0) dp.update_lz(i, 2, 2, {ind, 2}, Constant<1>(), pre16);
      }
    }
    if (phase + 1 < phase_total) {
      for (const auto c : candidate) pre[c] = 0;
      for (const auto cmd : dp.commands()) {
        if (cmd.type != pre16) continue;
        pre[candidate[cmd.val()]] += 1;
      }
      const size_t next_k = num_candidates[phase + 1];
      std::partial_sort(
        candidate.begin(), candidate.begin() + next_k, candidate.end(),
        [&](const uint16_t a, const uint16_t b) { return pre[a] > pre[b]; });
      for (size_t i = next_k; i < candidate.size(); ++i) pre[candidate[i]] = -1;
      candidate.resize(next_k);
    } else {
      using namespace data_type;
      writer ret;
      for (size_t i = 0; i < 64; ++i) ret.write<d16>(0);
      size_t adr = 0;
      for (const auto cmd : dp.commands()) {
        switch (cmd.type) {
        case uncomp: ret.write<d8, d8n>(cmd.len, {cmd.len, &input[adr]}); break;
        case rle: ret.write<d8, d8>(0x40 + cmd.len, input[adr]); break;
        case lz:
        case lzl: ret.write<d8, d16>(0x80 + cmd.len, cmd.lz_ofs); break;
        case pre16: ret.write<d8>(0xc0 + cmd.val()); break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      for (size_t i = 0; i < 64; ++i) write16(ret.out, 2 * i, candidate[i]);
      assert(dp.total_cost() + 0x80 == ret.size());
      assert(adr == input.size());
      ret.write<d8>(0);
      return ret.out;
    }
  }

  throw std::logic_error("phase_total == 0");
}

} // namespace sfc_comp
