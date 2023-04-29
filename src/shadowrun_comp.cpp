#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> shadowrun_comp_core(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  enum Tag {
    uncomp, lz
  };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      return li == rhs.li;
    }
    Tag tag;
    size_t li;
  };

  static constexpr size_t uncomp_max_len_bits = 15;
  static constexpr size_t lz_max_len_bits = 15;

  std::array<size_t, uncomp_max_len_bits + 1> uncomp_max_lens = {};
  for (size_t i = 0; i < uncomp_max_lens.size(); ++i) uncomp_max_lens[i] = (2 << i) - 1;
  std::array<size_t, lz_max_len_bits + 1> lz_max_lens;
  for (size_t i = 0; i < lz_max_lens.size(); ++i) lz_max_lens[i] = (2 << i) + 1;

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost = dp[i].cost;
    u_helper.update(i, cost);
    for (size_t k = 0; k < uncomp_max_lens.size(); ++k) {
      const size_t min_len = (k == 0) ? uncomp_max_lens[0] : uncomp_max_lens[k - 1] + 1;
      const size_t max_len = uncomp_max_lens[k];
      const auto res_u = u_helper.find(i + 1, min_len, max_len);
      dp.update_u(i + 1, res_u.len, {uncomp, k}, res_u.cost + 1 + (1 + 2 * k));
    }
    const auto res_lz = lz_helper.find_best(i, input.size());
    for (size_t k = 0; k < lz_max_lens.size(); ++k) {
      const size_t min_len = (k == 0) ? lz_max_lens[0] : lz_max_lens[k - 1] + 1;
      if (res_lz.len < min_len) break;
      const size_t max_len = lz_max_lens[k];
      const size_t base_cost = cost + 1 + (1 + 2 * k) + ilog2(2 * i + 1);
      dp.update_lz(i, min_len, max_len, res_lz, Constant<0>(), {lz, k}, base_cost);
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(4);
  writer_b8_h flags;

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      flags.write<b1>(false);
      for (size_t s = std::bit_floor(cmd.len) >> 1; s > 0; s >>= 1) {
        flags.write<b1, b1>(false, (cmd.len & s) != 0);
      }
      flags.write<b1>(true);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      assert(adr > 0);
      flags.write<b1>(true);
      flags.write<bnh>({ilog2(2 * adr + 1), cmd.lz_ofs});
      const size_t l = (cmd.len - lz_max_lens[0]) + 1;
      for (size_t s = std::bit_floor(l) >> 1; s > 0; s >>= 1) {
        flags.write<b1, b1>(false, (l & s) != 0);
      }
      flags.write<b1>(true);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 2);
  assert(adr == input.size());
  assert(dp.total_cost() + 32 == flags.bit_length() + 8 * ret.size());
  std::copy(flags.out.begin(), flags.out.end(), std::back_inserter(ret.out));

  return ret.out;
}

} // namespace

std::vector<uint8_t> shadowrun_comp(std::span<const uint8_t> input) {
  return shadowrun_comp_core(input);
}

std::vector<uint8_t> battletech_comp(std::span<const uint8_t> input) {
  auto ret = shadowrun_comp_core(input);
  ret.resize(ret.size() + 2);
  for (size_t i = ret.size() - 1; i >= 4; --i) ret[i] = ret[i - 2];
  write16(ret, 2, ret.size() - 4);
  return ret;
}

} // namespace sfc_comp
