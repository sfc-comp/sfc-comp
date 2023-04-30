#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> lennus_2_comp(std::span<const uint8_t> in) {
  enum CompType { uncomp, lz2, lz2_0, lz };

  static constexpr size_t pad = 0x800;

  // cf. $CB:8045
  static constexpr auto lz_lens= std::to_array<size_t>({
    3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 14, 16, 24, 32, 64, 128
  });

  std::array<size_t, lz_lens.back() + 1> len_code = {};
  for (size_t i = 0; i < lz_lens.size(); ++i) len_code[lz_lens[i]] = i;

  std::vector<uint8_t> input(in.size() + pad, 0);
  std::copy(in.begin(), in.end(), input.begin() + pad);

  std::vector<uint8_t> best;

  for (size_t comp_type = 0x5059; comp_type <= 0x5159; comp_type += 0x100) {
    lz_helper lz_helper(input);
    sssp_solver<CompType> dp(input.size(), pad);
    for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);

    for (size_t i = pad; i < input.size(); ++i) {
      dp.update(i, 1, 1, Constant<9>(), uncomp);
      if (comp_type == 0x5059) {
        const auto res_lz = lz_helper.find_best(i, 0x1000);
        dp.update_lz_table(i, lz_lens, res_lz, Constant<17>(), lz);
      } else {
        const auto res_lzs = lz_helper.find_best(i, 0x7f);
        dp.update_lz(i, 2, 2, res_lzs, Constant<9>(), lz2);
        if (i + 1 < input.size() && read16(input, i - 0x800) == read16(input, i)) {
          dp.update(i, 2, 2, Constant<9>(), lz2_0);
        }
        const auto res_lz = lz_helper.find_best(i, 0x800);
        dp.update_lz_table(i, lz_lens, res_lz, Constant<17>(), lz);
      }
      lz_helper.add_element(i);
    }

    using namespace data_type;
    writer_b8_h ret(8);
    size_t adr = pad;
    for (const auto& cmd : dp.commands(pad)) {
      const size_t d = adr - cmd.lz_ofs;
      switch (cmd.type) {
      case uncomp: ret.write<b1, d8>(false, input[adr]); break;
      case lz2_0: ret.write<b1, d8>(true, 0x80 | 0); break;
      case lz2: ret.write<b1, d8>(true, 0x80 | d); break;
      case lz: {
        if (comp_type == 0x5059) {
          ret.write<b1, d16>(true, (d & 0x0fff) << 4 | len_code[cmd.len]); break;
        } else {
          ret.write<b1, d16b>(true, (d & 0x07ff) << 4 | len_code[cmd.len]); break;
        }
      }
      default: assert(0);
      }
      adr += cmd.len;
    }
    write16(ret.out, 0, comp_type);
    write24(ret.out, 2, in.size());
    write24(ret.out, 5, ret.size() - 8);
    assert(adr == input.size());
    assert(dp.total_cost() + 8 * 8 == ret.bit_length());
    if (best.empty() || ret.size() < best.size()) best = std::move(ret.out);
  }
  return best;
}

} // namespace sfc_comp
