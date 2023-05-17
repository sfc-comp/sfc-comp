#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> sd_gundam_gnext_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0x8000);

  enum tag { uncomp, lz };

  const size_t pad = 0x22;
  std::vector<uint8_t> input(in.size() + pad);
  std::ranges::copy(in, input.begin() + pad);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0_2 = dp.c<0, 2>(0x21);

  for (size_t i = input.size(); i-- > pad; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    dp.update(i, 3, 0x21, lz_helper.find(i, 0x1000, 3), c0_2, 17, lz);
    c0_2.update(i);
  }

  using namespace data_type;
  writer_b8_l ret(2);
  size_t adr = pad;
  for (const auto& cmd : dp.optimal_path(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<b1, d8>(true, input[adr]); break;
    case lz: {
      assert(cmd.len & 1);
      uint16_t d = (cmd.lz_ofs() - pad) & 0x0fff;
      ret.write<b1, d8, d8>(false, d & 0x00ff, (d & 0x0f00) >> 8 | (cmd.len - 3) << 3);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.optimal_cost(pad) + 2 * 8 == ret.bit_length());
  assert(adr - pad == in.size());

  if (ret.size() < in.size() + 2) {
    write16(ret.out, 0, (in.size() - 1) | 0x8000);
  } else {
    write16(ret.out, 0, (in.size() - 1) | 0x0000);
    ret.out.resize(in.size() + 2);
    std::ranges::copy(in, ret.out.begin() + 2);
  }
  return ret.out;
}

} // namespace sfc_comp
