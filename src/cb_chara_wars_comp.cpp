#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> cb_chara_wars_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);

  enum CompType {
    uncomp, lz, lzl
  };

  lz_helper lz_helper(input);
  std::vector<encode::lz_data> lz_memo(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    lz_memo[i] = lz_helper.find_best(i, 0x1000);
    lz_helper.add_element(i);
  }

  // [TODO] Find a more clever way.
  const std::vector<size_t> lz_min_len_table = {2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 127};

  std::vector<uint8_t> best;
  for (size_t long_len = 0; long_len < 2; long_len++) {
    for (const size_t lz_min_len : lz_min_len_table) {
      sssp_solver<CompType> dp(input.size());

      for (size_t i = 0; i < input.size(); ++i) {
        dp.update(i, 1, 1, Constant<9>(), uncomp);
        auto res_lz = lz_memo[i];
        if (long_len) {
          dp.update_lz(i, lz_min_len, lz_min_len + 0x0e, res_lz, Constant<17>(), lz);
          dp.update_lz(i, lz_min_len + 0x0f, lz_min_len + 0x0f + 0xff, res_lz, Constant<25>(), lzl);
        } else {
          dp.update_lz(i, lz_min_len, lz_min_len + 0x0f, res_lz, Constant<17>(), lz);
        }
        lz_helper.add_element(i);
      }

      using namespace data_type;
      writer_b ret;
      size_t adr = 0;
      ret.write<d8, d16>((long_len) << 7 | lz_min_len, 0);
      for (const auto cmd : dp.commands()) {
        const size_t d = adr - cmd.lz_ofs;
        switch (cmd.type) {
        case uncomp: ret.write<b1l, d8>(true, input[adr]); break;
        case lz: ret.write<b1l, d16>(false, (cmd.len - lz_min_len) << 12 | (d - 1)); break;
        case lzl: ret.write<b1l, d16, d8>(false, 0xf000 | (d - 1), (cmd.len - lz_min_len - 0x0f)); break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert((dp.total_cost() + 7) / 8 + 3 == ret.size());
      assert(adr == input.size());
      write16(ret.out, 1, input.size());
      if (best.empty() || ret.size() < best.size()) {
        best = std::move(ret.out);
      }
    }
  }

  return best;
}

} // namespace sfc_comp
