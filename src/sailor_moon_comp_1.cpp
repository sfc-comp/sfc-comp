#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sailor_moon_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum CompType {
    uncomp, lzs, lzls, lzll
  };
  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lzs = lz_helper.find_best(i, 0x100);
    dp.update_lz(i, 2, 5, res_lzs, Constant<12>(), lzs);
    auto res_lzl = lz_helper.find_best(i, 0x2000);
    dp.update_lz(i, 3, 9, res_lzl, Constant<18>(), lzls);
    dp.update_lz(i, 10, 256, res_lzl, Constant<26>(), lzll);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  ret.write<d16>(0);

  size_t curr16 = 0, bit16_pos = 0;
  auto write_b16 = [&](size_t input_bits, size_t v) {
    while (input_bits > 0) {
      --input_bits;
      if ((v >> input_bits) & 1) {
        ret.out[curr16 + (bit16_pos >> 3)] |= 1 << (bit16_pos & 7);
      }
      ++bit16_pos;
      if (bit16_pos == 16) {
        bit16_pos = 0;
        curr16 = ret.size();
        ret.write<d16>(0);
      }
    }
  };

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;

    switch (cmd.type) {
    case uncomp: {
      write_b16(1, 1);
      ret.write<d8>(input[adr]);
    } break;
    case lzs: {
      write_b16(4, cmd.len - 2);
      ret.write<d8>(0x100 - d);
    } break;
    case lzls: case lzll: {
      write_b16(2, 1);
      ret.write<d8>((0x2000 - d) & 0x00ff);
      if (cmd.type == lzls) {
        ret.write<d8>((cmd.len - 2) | ((0x2000 - d) & 0x1f00) >> 5);
      } else {
        ret.write<d8>(((0x2000 - d) & 0x1f00) >> 5);
        ret.write<d8>(cmd.len - 1);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write_b16(2, 1);
  ret.write<d16, d8>(0xf000, 0);

  size_t cost = (dp.total_cost() + 2 + 7) / 8 + 3;
  assert(cost <= ret.size() && ret.size() <= cost + 1 + (bit16_pos == 0 ? 1 : 0));
  return ret.out;
}

} // namespace sfc_comp
