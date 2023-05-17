#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> cb_chara_wars_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);

  enum tag { uncomp, lz, lzl };

  // [TODO] Find a clever way.
  static constexpr auto min_lens = std::to_array<size_t>({
    2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 127
  });

  lz_helper lz_helper(input);
  std::vector<encode::lz_data> lz_memo(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    lz_memo[i] = lz_helper.find(i, 0x1000, min_lens.front());
    lz_helper.add_element(i);
  }

  std::vector<uint8_t> best;
  for (size_t t = 0; t < 2; t++) {
    bool long_len = t > 0;
    for (const size_t min_len : min_lens) {
      solver<tag> dp(input.size());
      auto c0 = dp.c<0>(min_len + 0x0f + (long_len ? 0xff : 0));

      for (size_t i = input.size(); i-- > 0; ) {
        dp.update(i, 1, 9, uncomp);
        const auto res_lz = lz_memo[i];
        if (long_len) {
          dp.update(i, min_len, min_len + 0x0e, res_lz, c0, 17, lz);
          dp.update(i, min_len + 0x0f, min_len + 0x0f + 0xff, res_lz, c0, 25, lzl);
        } else {
          dp.update(i, min_len, min_len + 0x0f, res_lz, c0, 17, lz);
        }
        c0.update(i);
      }

      using namespace data_type;
      writer_b8_l ret; ret.write<d8, d16>(t << 7 | min_len, 0);
      size_t adr = 0;
      for (const auto& cmd : dp.optimal_path()) {
        const size_t d = adr - cmd.lz_ofs();
        switch (cmd.type) {
        case uncomp: ret.write<b1, d8>(true, input[adr]); break;
        case lz: ret.write<b1, d16>(false, (cmd.len - min_len) << 12 | (d - 1)); break;
        case lzl: ret.write<b1, d16, d8>(false, 0xf000 | (d - 1), (cmd.len - min_len - 0x0f)); break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());
      assert(dp.optimal_cost() + 3 * 8 == ret.bit_length());
      write16(ret.out, 1, input.size());
      if (best.empty() || ret.size() < best.size()) best = std::move(ret.out);
    }
  }

  return best;
}

} // namespace sfc_comp
