#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> soccer_kid_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);
  enum tag { uncomp, uncompl, lz2, lz3, lz4, lz };

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    u_helper.update(i, dp[i].cost);
    const auto u1 = u_helper.find(i + 1, 1, 8);
    dp.update_u(i + 1, u1.len, uncomp, u1.cost + 5);
    const auto u2 = u_helper.find(i + 1, 9, 9 + 0xff);
    dp.update_u(i + 1, u2.len, uncompl, u2.cost + 11);

    const auto res_lz2 = lz_helper.find(i, 0xff, 2);
    dp.update_lz(i, 2, 2, res_lz2, Constant<10>(), lz2);
    const auto res_lz3 = lz_helper.find(i, 0x1ff, 3);
    dp.update_lz(i, 3, 3, res_lz3, Constant<12>(), lz3);
    const auto res_lz4 = lz_helper.find(i, 0x3ff, 4);
    dp.update_lz(i, 4, 4, res_lz4, Constant<13>(), lz4);
    const auto res_lz = lz_helper.find(i, 0xfff, 5);
    dp.update_lz(i, 5, 256, res_lz, Constant<23>(), lz);

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_l ret(4); ret.write<b1>(false);

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp:
    case uncompl: {
      if (cmd.type == uncomp) ret.write<bnh, bnh>({2, 0}, {3, cmd.len - 1});
      else ret.write<bnh, bnh>({3, 7}, {8, cmd.len - 9});
      ret.write<b8hn>({cmd.len, &input[adr]});
    } break;
    case lz2: ret.write<bnh, bnh>({2, 1}, {8, d}); break;
    case lz3: ret.write<bnh, bnh>({3, 4}, {9, d}); break;
    case lz4: ret.write<bnh, bnh>({3, 5}, {10, d}); break;
    case lz: ret.write<bnh, bnh, bnh>({3, 6}, {8, cmd.len - 1}, {12, d}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 1 + 32 == ret.bit_length());

  const size_t s = std::min<size_t>(8, ret.size());
  uint64_t v = 1;
  for (size_t i = s - 1; i >= 4; --i) v = v << 8 | ret[i];
  v >>= 1;
  for (size_t i = 4; i < s; ++i) ret[i] = v & 0xff, v >>= 8;

  uint32_t checksum = 0;
  for (size_t i = 4; i + 3 < ret.size(); i += 4) checksum ^= read32(ret.out, i);
  for (size_t r = ret.size() % 4, j = 0; j < r; ++j) checksum ^= ret[ret.size() - r + j] << (8 * j);

  std::reverse(ret.out.begin() + 4, ret.out.end());
  write32b(ret.out, 0, ret.size() - 4);
  ret.write<d32b, d32b>(checksum, input.size());

  return ret.out;
}

} // namespace sfc_comp
