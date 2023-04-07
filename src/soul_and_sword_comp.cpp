#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> soul_and_sword_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 2, 0x10000);

  enum Tag {
    uncomp, lz
  };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag == uncomp) {
        return len_no;
      } else {
        return len_no == rhs.len_no && ofs_no == rhs.ofs_no; // ?
      }
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  static constexpr size_t uncomp_len_tab[][3] = {
    {0x0001, 0x0000,  3},
    {0x0004, 0x0006,  4},
    {0x0005, (((1 << 5) - 3) << 0),  6},
    {0x0006, (((1 << 6) - 3) << 1),  8},
    {0x0008, (((1 << 7) - 3) << 2),  10},
    {0x000c, (((1 << 8) - 3) << 3),  12},
    {0x0014, (((1 << 9) - 3) << 4),  14},
    {0x0024, (((1 << 10) - 3) << 5),  16},
    {0x0044, (((1 << 11) - 3) << 6),  18},
    {0x0084, (((1 << 12) - 3) << 7),  20},
    {0x0101, 0,  0}
  };

  static constexpr size_t lz_len_tab[][3] = {
    {0x0002, 0x0000,  3},
    {0x0009, 0x000e,  4},
    {0x000a, (((1 << 6) - 3) << 0),  6},
    {0x000b, (((1 << 7) - 3) << 1),  8},
    {0x000d, (((1 << 8) - 3) << 2),  10},
    {0x0011, (((1 << 9) - 3) << 3),  12},
    {0x0019, (((1 << 10) - 3) << 4),  14},
    {0x0029, (((1 << 11) - 3) << 5),  16},
    {0x0049, (((1 << 12) - 3) << 6),  18},
    {0x0089, (((1 << 13) - 3) << 7),  20},
    {0x0100, 0, 0},
  };

  for (size_t i = 0; i < 2; ++i) lz_helper.add_element(i);
  for (size_t i = 0; i < 2; ++i) dp[i + 1].cost = 0;

  for (size_t i = 2; i < input.size(); ++i) {
    const auto cost = dp[i].cost;

    for (size_t k = 0; k < 10; ++k) {
      const size_t min_len = uncomp_len_tab[k][0];
      if (i + min_len > input.size()) break;
      const size_t max_len = uncomp_len_tab[k + 1][0] - 1;
      const size_t len_bitsize = uncomp_len_tab[k][2];
      const size_t total_cost = cost + len_bitsize;
      dp.update(i, min_len, max_len, Linear<8, 0>(), {uncomp, 0, k}, total_cost);
    }

    const size_t ofs_bitsize = std::min<size_t>(12, 1 + ilog2(i));
    auto res_lz = lz_helper.find_best(i, 0x0fff);

    for (size_t k = 0; k < 10; ++k) {
      const size_t min_len = lz_len_tab[k][0];
      if (min_len > res_lz.len) break;
      const size_t max_len = lz_len_tab[k + 1][0] - 1;
      const size_t len_bitsize = lz_len_tab[k][2];
      const size_t total_cost = cost + 1 + ofs_bitsize + len_bitsize;
      dp.update_lz(i, min_len, max_len, res_lz,
                   Constant<0>(), {lz, ofs_bitsize, k}, total_cost);
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  writer raws;
  ret.write<d16, d16>(0, 0);
  raws.write<d8, d8>(input[0], input[1]);

  size_t adr = 2;
  for (const auto cmd : dp.commands(2)) {
    switch (cmd.type.tag) {
    case uncomp: {
      const size_t li = cmd.type.len_no;
      ret.write<b8hn_h>({uncomp_len_tab[li][2], uncomp_len_tab[li][1] + (cmd.len - uncomp_len_tab[li][0])});
      raws.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const size_t li = cmd.type.len_no;
      ret.write<b1h>(true);
      assert(d < (size_t(1) << cmd.type.ofs_no));
      ret.write<b8hn_h>({cmd.type.ofs_no, d});
      ret.write<b8hn_h>({lz_len_tab[li][2], lz_len_tab[li][1] + (cmd.len - lz_len_tab[li][0])});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 4);
  std::copy(raws.out.begin(), raws.out.end(), std::back_inserter(ret.out));

  assert(ret.size() == 4 + 2 + (dp.total_cost() + 7) / 8);

  return ret.out;
}

} // namespace sfc_comp
