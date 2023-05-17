#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> soccer_kid_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);
  enum tag { uncomp, uncompl, lz2, lz3, lz4, lz };

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(256);
  auto c8 = dp.c<8>(9 + 0xff);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 8, c8, 5, uncomp);
    dp.update(i, 9, 9 + 0xff, c8, 11, uncompl);

    dp.update(i, 2, 2, lz_helper.find(i, 0xff, 2), c0, 10, lz2);
    dp.update(i, 3, 3, lz_helper.find(i, 0x1ff, 3), c0, 12, lz3);
    dp.update(i, 4, 4, lz_helper.find(i, 0x3ff, 4), c0, 13, lz4);
    dp.update(i, 5, 256, lz_helper.find(i, 0xfff, 5), c0, 23, lz);

    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b8_l ret(4); ret.write<b1>(false);

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp:
    case uncompl: {
      if (cmd.type == uncomp) ret.write<bnh, bnh>({2, 0b00}, {3, cmd.len - 1});
      else ret.write<bnh, bnh>({3, 0b111}, {8, cmd.len - 9});
      ret.write<b8hn>({cmd.len, &input[adr]});
    } break;
    case lz2: ret.write<bnh, bnh>({2, 0b01}, {8, d}); break;
    case lz3: ret.write<bnh, bnh>({3, 0b100}, {9, d}); break;
    case lz4: ret.write<bnh, bnh>({3, 0b101}, {10, d}); break;
    case lz: ret.write<bnh, bnh, bnh>({3, 0b110}, {8, cmd.len - 1}, {12, d}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 1 + 32 == ret.bit_length());

  const size_t s = std::min<size_t>(8, ret.size());
  uint64_t v = 1;
  for (size_t i = s - 1; i >= 4; --i) v = v << 8 | ret[i];
  v >>= 1;
  for (size_t i = 4; i < s; ++i) ret[i] = v & 0xff, v >>= 8;

  uint32_t checksum = 0;
  for (size_t i = 4; i + 3 < ret.size(); i += 4) checksum ^= read32(ret.out, i);
  for (size_t r = ret.size() % 4, j = 0; j < r; ++j) checksum ^= ret[ret.size() - r + j] << (8 * j);

  std::reverse(ret.out.begin() + 4, ret.out.end());
  ret.write<d32b, d32b>(checksum, input.size());
  write32b(ret.out, 0, ret.size() - 4);

  return ret.out;
}

} // namespace sfc_comp
