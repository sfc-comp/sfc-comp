#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> stargate_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);
  enum Tag {
    none, uncomp, lz
  };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == none) return true;
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

  if (input.size() > uncomp_max_lens.back()) {
    throw std::runtime_error("This algorithm cannot compress the given data.");
  }

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<CompType> dp0(input.size()), dp1(input.size());
  dp1.reset(0);

  for (size_t i = 0; i < input.size(); ++i) {
    const auto cost0 = dp0[i].cost;
    u_helper.update(i, cost0);

    dp1.update(i, 0, 0, Constant<1>(), {none, 0}, cost0);
    for (size_t k = 0; k < uncomp_max_lens.size(); ++k) {
      const size_t min_len = (k == 0) ? uncomp_max_lens[0] : uncomp_max_lens[k - 1] + 1;
      const size_t max_len = uncomp_max_lens[k];
      const auto res_u = u_helper.find(i + 1, min_len, max_len);
      dp1.update_u(i + 1, res_u.len, {uncomp, k}, res_u.cost + 1 + (1 + 2 * k));
    }
    const auto cost1 = dp1[i].cost;
    const auto res_lz = lz_helper.find_best(i, input.size());
    for (size_t k = 0; k < lz_max_lens.size(); ++k) {
      const size_t min_len = (k == 0) ? lz_max_lens[0] : lz_max_lens[k - 1] + 1;
      if (res_lz.len < min_len) break;
      const size_t max_len = lz_max_lens[k];
      const size_t base_cost = cost1 + (1 + 2 * k) + ilog2(2 * i + 1);
      dp0.update_lz(i, min_len, max_len, res_lz, Constant<0>(), {lz, k}, base_cost);
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret; ret.write<d16, d16>(0, 0);
  writer_b flags;

  using command_type = decltype(dp0)::vertex_type;
  const std::vector<command_type> commands = [&] {
    std::vector<command_type> ret;
    size_t adr = input.size();
    size_t curr = (dp0.total_cost() < dp1.total_cost()) ? 0 : 1;
    while (adr > 0) {
      command_type cmd;
      if (curr == 0) cmd = dp0[adr], curr = 1, assert(cmd.len > 0);
      else cmd = dp1[adr], curr = 0;
      adr -= cmd.len;
      ret.emplace_back(cmd);
    }
    assert(adr == 0 && curr == 0);
    std::reverse(ret.begin(), ret.end());
    return ret;
  }();

  size_t adr = 0;
  for (const auto cmd : commands) {
    switch (cmd.type.tag) {
    case none: {
      flags.write<b1h>(false);
    } break;
    case uncomp: {
      if (adr > 0) flags.write<b1h>(true);
      for (size_t s = (1 << ilog2(cmd.len)) >> 1; s > 0; s >>= 1) {
        flags.write<b1h, b1h>(true, (cmd.len & s) != 0);
      }
      flags.write<b1h>(false);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      assert(adr > 0);
      flags.write<b8hn_h>({ilog2(2 * adr + 1), cmd.lz_ofs});
      const size_t l = (cmd.len - lz_max_lens[0]) + 1;
      for (size_t s = (1 << ilog2(l)) >> 1; s > 0; s >>= 1) {
        flags.write<b1h, b1h>(true, (l & s) != 0);
      }
      flags.write<b1h>(false);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 2);
  assert(adr == input.size());
  assert((std::min(dp0.total_cost(), dp1.total_cost()) - 1) + 32 == flags.bit_length() + 8 * ret.size());
  std::copy(flags.out.begin(), flags.out.end(), std::back_inserter(ret.out));

  return ret.out;
}

} // namespace sfc_comp
