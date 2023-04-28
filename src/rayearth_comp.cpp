#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> rayearth_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum CompType {
    uncomp, uncompl,
    lz, lzl,
    lz4, lz4l
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 0x01, 0x3f, Linear<1, 1>(), uncomp);
    dp.update(i, 0x40, 0xff, Linear<1, 2>(), uncompl);
    if (auto res_lz4 = lz_helper.find_best(i, 0x03fe); res_lz4.len >= 3) {
      if (i - res_lz4.ofs >= 0x02ff) {
        dp.update_lz(i, 0x0003, 0x0041, res_lz4, Constant<3>(), lz4);
        dp.update_lz(i, 0x0042, 0x00ff, res_lz4, Constant<4>(), lz4l);
      }
      if (auto res_lz = lz_helper.find_best(i, 0x02fe); res_lz.len >= 3) {
        dp.update_lz(i, 0x0003, 0x0041, res_lz, Constant<2>(), lz);
        dp.update_lz(i, 0x0042, 0x00ff, res_lz, Constant<3>(), lzl);
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
    case lz: case lzl: {
      if (cmd.type == lz) ret.write<d8>(0x40 + ((d & 0x300) >> 2 | (cmd.len - 3)));
      else ret.write<d8, d8>(0x40 + ((d & 0x300) >> 2 | 0x3f), cmd.len - 0x42);
      ret.write<d8>(d & 0xff);
    } break;
    case lz4: case lz4l: {
      if (cmd.type == lz4) ret.write<d8>(0xc0 | (cmd.len - 3));
      else ret.write<d8, d8>(0xff, cmd.len - 0x42);
      ret.write<d8, d8>(0xff, d - 0x2ff);
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
