#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> rayearth_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum CompType {
    uncomp, uncompl,
    lz1, lz1l,
    lz2, lz2l,
    lz3, lz3l,
    lz4, lz4l
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 0x01, 0x3f, Linear<1, 1>(), uncomp);
    dp.update(i, 0x40, 0xff, Linear<1, 2>(), uncompl);

    auto res_lz4 = lz_helper.find_best(i, 0x03fe);
    if (res_lz4.len >= 3) {
      if (i - res_lz4.ofs >= 0x02ff) {
        dp.update_lz(i, 0x0003, 0x0041, res_lz4, Constant<3>(), lz4);
        dp.update_lz(i, 0x0042, 0x00ff, res_lz4, Constant<4>(), lz4l);
      }
      auto res_lz = lz_helper.find_best(i, 0x02fe);
      if (res_lz.len >= 3) {
        const size_t d = i - res_lz.ofs;
        if (d >= 0x0200) {
          dp.update_lz(i, 0x0003, 0x0041, res_lz, Constant<2>(), lz3);
          dp.update_lz(i, 0x0042, 0x00ff, res_lz, Constant<3>(), lz3l);
        } else if (d >= 0x0100) {
          dp.update_lz(i, 0x0003, 0x0041, res_lz, Constant<2>(), lz2);
          dp.update_lz(i, 0x0042, 0x00ff, res_lz, Constant<3>(), lz2l);
        } else {
          dp.update_lz(i, 0x0003, 0x0041, res_lz, Constant<2>(), lz1);
          dp.update_lz(i, 0x0042, 0x00ff, res_lz, Constant<3>(), lz1l);
        }
      }
    }

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret;
  ret.write<d8, d16>(0, 0);

  size_t adr = 0;
  for (const auto cmd : dp.commands(0)) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: case uncompl: {
      if (cmd.type == uncomp) ret.write<d8>(cmd.len - 1);
      else ret.write<d8, d8>(0x3f, cmd.len - 0x40);
      ret.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz1: case lz1l: {
      if (cmd.type == lz1) ret.write<d8>(0x40 | (cmd.len - 3));
      else ret.write<d8, d8>(0x7f, cmd.len - 0x42);
      ret.write<d8>(d);
    } break;
    case lz2: case lz2l: {
      if (cmd.type == lz2) ret.write<d8>(0x80 | (cmd.len - 3));
      else ret.write<d8, d8>(0xbf, cmd.len - 0x42);
      ret.write<d8>(d - 0x100);
    } break;
    case lz3: case lz3l:
    case lz4: case lz4l: {
      if (cmd.type == lz3 || cmd.type == lz4) ret.write<d8>(0xc0 | (cmd.len - 3));
      else ret.write<d8, d8>(0xff, cmd.len - 0x42);
      if (cmd.type == lz3 || cmd.type == lz3l) ret.write<d8>(d - 0x200);
      else ret.write<d8, d8>(0xff, d - 0x2ff);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.out[0] = 0;
  write16(ret.out, 1, input.size());
  assert(ret.size() == 3 + dp.total_cost());
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
