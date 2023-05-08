#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> stargate_comp_core(std::span<const uint8_t> input, const bool terminator_b, const bool uncomp_b) {
  check_size(input.size(), 1, 0xffff);

  enum method { none, uncomp, lz };
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

  if (input.size() > uncomp_lens.back().max) {
    throw std::runtime_error("This algorithm may not be able to compress the given data.");
  }

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp0(input.size()), dp1(input.size(), -1);

  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost0 = dp0[i].cost;
    u_helper.update(i, cost0);

    dp1.update(i, 0, 0, Constant<1>(), {none, 0}, cost0);
    for (size_t k = 0; k < uncomp_lens.size(); ++k) {
      const auto res_u = u_helper.find(i + 1, uncomp_lens[k].min, uncomp_lens[k].max);
      if (res_u.len == u_helper.nlen) continue;
      dp1.update_u(i + 1, res_u.len, {uncomp, k}, res_u.cost + 1 + uncomp_lens[k].bitlen);
    }
    dp0.update_lz_matrix(i, lz_offsets, lz_lens,
      [&](size_t oi) { return lz_helper.find_best(i, lz_offsets[oi].max); },
      [&](size_t, size_t li) -> tag { return {lz, li}; },
      ilog2(2 * i + 1), dp1[i].cost
    );
    lz_helper.add_element(i);
  }

  const auto [commands, min_cost] = [&]{
    using command_type = decltype(dp0)::vertex_type;
    std::vector<command_type> ret;
    const auto min_cost = std::min(dp0.total_cost(), dp1.total_cost());
    ptrdiff_t adr = input.size();
    size_t curr = (dp0.total_cost() == min_cost) ? 0 : 1;
    while (adr > 0 || (adr == 0 && curr > 0)) {
      command_type cmd;
      if (curr == 0) cmd = dp0[adr], curr = 1, assert(cmd.len > 0);
      else cmd = dp1[adr], curr = 0;
      adr -= cmd.len;
      ret.emplace_back(cmd);
    }
    assert(adr == 0 && curr == 0);
    std::reverse(ret.begin(), ret.end());
    return std::make_pair(std::move(ret), min_cost);
  }();

  using namespace data_type;
  writer ret(4);
  writer_b8_h flags;

  const auto write_len = [&](size_t l) {
    size_t s = std::bit_floor(l) >> 1;
    for (; s > 0; s >>= 1) flags.write<b1, b1>(!terminator_b, (l & s) != 0);
    flags.write<b1>(terminator_b);
  };

  size_t adr = 0;
  for (const auto cmd : commands) {
    switch (cmd.type.tag) {
    case none: {
      flags.write<b1>(!uncomp_b);
    } break;
    case uncomp: {
      if (adr > 0) flags.write<b1>(uncomp_b);
      write_len(cmd.len);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      assert(adr > 0);
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
  assert((min_cost - 1) + 32 == flags.bit_length() + 8 * ret.size());
  ret.extend(flags);

  return ret.out;
}

} // namespace

std::vector<uint8_t> shadowrun_comp(std::span<const uint8_t> input) {
  return stargate_comp_core(input, true, false);
}

std::vector<uint8_t> stargate_comp(std::span<const uint8_t> input) {
  return stargate_comp_core(input, false, true);
}

} // namespace sfc_comp
