#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_4wd_the_baja_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);
  enum tag { uncomp, uncompl, lz2, lz3, lz4, lz };

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(5 + 0xff);
  auto c8 = dp.c<8>(9 + 0xff);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 8, c8, 5, uncomp);
    dp.update(i, 9, 9 + 0xff, c8, 11, uncompl);

    const auto res_lz = lz_helper.find(i, 0xff, 2);
    dp.update(i, 2, 2, lz_helper.find(i, 0xff, 2), c0, 10, lz2);
    dp.update(i, 3, 3, lz_helper.find(i, 0x1ff, 3), c0, 12, lz3);
    dp.update(i, 4, 4, lz_helper.find(i, 0x3ff, 4), c0, 13, lz4);
    dp.update(i, 5, 5 + 0xff, res_lz, c0, 19, lz);

    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b16_l ret(2); ret.write<b1>(false);

  const auto write8 = [&ret](size_t v) {
    if (ret.bit == 0 || ret.bit >= 8) ret.write<bnl>({8, v});
    else ret.write<bnh>({8, v});
  };

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp:
    case uncompl: {
      if (cmd.type == uncomp) ret.write<bnh, bnh>({2, 0b00}, {3, cmd.len - 1});
      else { ret.write<bnh>({3, 0b111}); write8(cmd.len - 9); }
      for (size_t i = 0; i < cmd.len; ++i) write8(input[adr + i]);
    } break;
    case lz2: { ret.write<bnh>({2, 0b01}); write8(d); } break;
    case lz3: { ret.write<bnh, bnh>({3, 0b100}, {9, d}); } break;
    case lz4: { ret.write<bnh, bnh>({3, 0b101}, {10, d}); } break;
    case lz:  { ret.write<bnh>({3, 0b110}); write8(cmd.len - 5); write8(d); } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 1 + 16 == ret.bit_length());

  if (ret.size() < 4) throw std::logic_error("Input is empty.");

  write16(ret.out, 2, 0x8000 | read16(ret.out, 2) >> 1);
  assert(ret.size() % 2 == 0);
  std::reverse(ret.out.begin() + 2, ret.out.end());
  for (size_t i = 2; i < ret.size(); i += 2) std::swap(ret[i], ret[i + 1]);
  ret.write<d16>(0);
  write16(ret.out, 0, ret.size() - 4);
  write16(ret.out, ret.size() - 2, input.size());

  return ret.out;
}

} // namespace sfc_comp
