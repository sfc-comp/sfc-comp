#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> tales_of_phantasia_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum CompType {
    uncomp, rle, rlel, lz
  };

  std::vector<uint8_t> best;
  for (size_t comp_type = 0x81; comp_type <= 0x83; comp_type += 2) {
    lz_helper lz_helper(input);
    sssp_solver<CompType> dp(input.size());

    size_t rlen = 0;
    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, Constant<9>(), uncomp);
      if (comp_type == 0x81) {
        auto res_lz = lz_helper.find_best(i, 0x0fff);
        dp.update_lz(i, 3, 0x12, res_lz, Constant<17>(), lz);
      } else {
        auto res_lz = lz_helper.find_best(i, 0x0fff);
        dp.update_lz(i, 3, 0x11, res_lz, Constant<17>(), lz);
        rlen = encode::run_length(input, i, rlen);
        dp.update(i, 4, 0x12, rlen, Constant<17>(), rle);
        dp.update(i, 0x13, 0x13 + 0xff, rlen, Constant<25>(), rlel);
      }
      lz_helper.add_element(i);
    }

    using namespace data_type;
    writer_b ret; ret.write<d8, d32, d32>(0, 0, 0);
    size_t adr = 0;
    for (const auto cmd : dp.commands()) {
      size_t d = adr - cmd.lz_ofs;
      switch (cmd.type) {
      case uncomp: ret.write<b1l, d8>(true, input[adr]); break;
      case lz: ret.write<b1l, d16>(false, d | (cmd.len - 3) << 12); break;
      case rle: ret.write<b1l, d8, d8>(false, input[adr], 0xf0 | ((cmd.len - 4) + 1)); break;
      case rlel: ret.write<b1l, d8, d8, d8>(false, cmd.len - 0x13, 0xf0, input[adr]); break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    ret.out[0] = comp_type;
    write32(ret.out, 1, ret.out.size() - 9);
    write32(ret.out, 5, input.size());
    assert((dp.total_cost() + 7) / 8 + 9 == ret.size());
    assert(adr == input.size());
    if (best.size() == 0 || ret.out.size() < best.size()) {
      best = std::move(ret.out);
    }
  }
  if (best.size() >= input.size() + 9) {
    best.resize(input.size() + 9);
    best[0] = 0x80;
    write32(best, 1, input.size());
    write32(best, 5, input.size());
    std::copy(input.begin(), input.end(), best.begin() + 9);
  }
  return best;
}

} // namespace sfc_comp
