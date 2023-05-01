#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> smash_tv_comp_core(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum method { uncomp, lz };
  using tag = tag_l<method>;

  static constexpr size_t uncomp_len_max_bits = 15;
  static constexpr size_t lz_len_max_bits = 15;
  static constexpr size_t lz_min_len = 3;

  static constexpr auto uncomp_lens = create_array<vrange, uncomp_len_max_bits + 1>([](size_t i) {
    return vrange(1 << i, (2 << i) - 1, 2 * i + 1, 0);
  });
  static constexpr auto lz_lens = create_array<vrange, lz_len_max_bits + 1>([](size_t i) {
    return vrange((1 << i) + lz_min_len - 1, (2 << i) + lz_min_len - 2, 2 * i + 1, 0);
  });
  const std::array<vrange, 1> lz_offsets = {vrange(1, input.size(), 0, 0)};

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost = dp[i].cost;
    u_helper.update(i, cost);
    for (size_t k = 0; k < uncomp_lens.size(); ++k) {
      const auto res_u = u_helper.find(i + 1, uncomp_lens[k].min, uncomp_lens[k].max);
      dp.update_u(i + 1, res_u.len, {uncomp, k}, res_u.cost + 1 + uncomp_lens[k].bitlen);
    }
    dp.update_lz_matrix(i, lz_offsets, lz_lens,
      [&](size_t oi) { return lz_helper.find_best(i, lz_offsets[oi].max); },
      [&](size_t, size_t li) -> tag { return {lz, li}; },
      1 + ilog2(2 * i + 1)
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(4);
  writer_b8_h flags;

  const auto write_len = [&](size_t l) {
    size_t s = std::bit_floor(l) >> 1;
    for (; s > 0; s >>= 1) flags.write<b1, b1>(false, (l & s) != 0);
    flags.write<b1>(true);
  };

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      flags.write<b1>(false);
      write_len(cmd.len);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      assert(adr > 0);
      flags.write<b1>(true);
      flags.write<bnh>({ilog2(2 * adr + 1), cmd.lz_ofs});
      write_len(cmd.len - (lz_min_len - 1));
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

std::vector<uint8_t> smash_tv_comp(std::span<const uint8_t> input) {
  return smash_tv_comp_core(input);
}

std::vector<uint8_t> battletech_comp(std::span<const uint8_t> input) {
  auto ret = smash_tv_comp_core(input);
  ret.resize(ret.size() + 2);
  for (size_t i = ret.size() - 1; i >= 4; --i) ret[i] = ret[i - 2];
  write16(ret, 2, ret.size() - 4);
  return ret;
}

} // namespace sfc_comp
