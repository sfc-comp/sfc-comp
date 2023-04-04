#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> chrono_trigger_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  enum CompType {
    uncomp, lz
  };

  std::vector<uint8_t> best;

  for (size_t ty = 0; ty < 2; ++ty) {
    const size_t lz_max_len = (0x0010 << ty) - 1 + 3;
    const size_t lz_max_ofs = (0x1000 >> ty) - 1;

    lz_helper lz_helper(input);
    sssp_solver<CompType> dp(input.size());

    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, Constant<9>(), uncomp);
      auto res_lz = lz_helper.find_best(i, lz_max_ofs);
      dp.update_lz(i, 3, lz_max_len, res_lz, Constant<17>(), lz);
      lz_helper.add_element(i);
    }

    using namespace data_type;
    writer_b ret;
    size_t adr = 0;
    ret.write<d16>(0);
    for (const auto cmd : dp.commands()) {
      switch (cmd.type) {
      case uncomp: ret.write<b1l, d8>(false, input[adr]); break;
      case lz: {
        uint16_t d = adr - cmd.lz_ofs;
        ret.write<b1l, d16>(true, d | (cmd.len - 3) << (12 - ty));
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    assert((dp.total_cost() + 7) / 8 + 2 == ret.size());
    assert(adr == input.size());

    size_t comp_type_bit = (ty == 0) ? 0x00 : 0x40;
    if (ret.bit == 0) {
      write16(ret.out, 0, ret.size() - 2);
    } else {
      write16(ret.out, 0, ret.bits_pos - 2);
      ret.write<d24>(0);
      size_t len = (8 - ret.bit) + popcount32(ret.out[ret.bits_pos]) + 1;
      assert(ret.size() == ret.bits_pos + 3 + len);
      for (size_t i = 0, o = ret.bits_pos + len - 1; i < len; ++i) {
        ret.out[o + 3 - i] = ret.out[o - i];
      }
      ret.out[ret.bits_pos] = (8 - ret.bit) | comp_type_bit;
      write16(ret.out, ret.bits_pos + 1, ret.size());
    }
    ret.write<d8>(comp_type_bit);
    if (best.size() == 0 || ret.size() < best.size()) {
      best = std::move(ret.out);
    }
  }
  return best;
}

} // namespace sfc_comp
