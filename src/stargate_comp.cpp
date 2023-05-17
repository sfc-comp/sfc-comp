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

  if (input.size() > uncomp_lens.back().max) {
    throw std::runtime_error("This algorithm may not be able to compress the given data.");
  }

  lz_helper lz_helper(input, true);
  solver<tag> dp0(input.size()), dp1(input.size());
  auto c0_0 = dp0.c<0>(lz_lens.back().max);
  auto c8_1 = dp1.c<8>(uncomp_lens.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    const auto res_lz = lz_helper.find(i, input.size(), lz_lens.front().min);
    dp1.update(i, lz_lens, res_lz, c0_0, ilog2(2 * i + 1),
      [&](size_t li) -> tag { return {lz, li}; }
    );
    dp0.update(i, uncomp_lens, c8_1, 1, [&](size_t li) -> tag {return {uncomp, li}; });
    dp0.update_c(i, 0, dp1[i].cost + 1, {none, 0});
    c0_0.update(i); c8_1.update(i);
  }

  using namespace data_type;
  writer ret(4);
  writer_b8_h flags;

  const auto write_len = [&](size_t l) {
    size_t s = std::bit_floor(l) >> 1;
    for (; s > 0; s >>= 1) flags.write<b1, b1>(!terminator_b, (l & s) != 0);
    flags.write<b1>(terminator_b);
  };

  size_t adr = 0;
  for (size_t curr = 0; adr < input.size(); curr ^= 1) {
    const auto& cmd = (curr == 0) ? dp0[adr] : dp1[adr];
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
      flags.write<bnh>({ilog2(2 * adr + 1), cmd.lz_ofs()});
      write_len(cmd.len - (lz_min_len - 1));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 2);
  assert(adr == input.size());
  assert((dp0.optimal_cost() - 1) + 32 == flags.bit_length() + 8 * ret.size());
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
