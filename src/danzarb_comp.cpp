#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> danzarb_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum tag { uncomp, lz };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(32);
  auto c8 = dp.c<8>(8);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 8, c8, 4, uncomp);
    if (i > 0) {
      dp.update(i, 1, 32, lz_helper.find(i, input.size(), 1), c0, 6 + ilog2(2 * i - 1), lz);
    }
    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(2);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    switch (cmd.type) {
    case uncomp: ret.write<b1, bnh, b8hn>(false, {3, cmd.len & 7}, {cmd.len, &input[adr]}); break;
    case lz: ret.write<b1, bnh, bnh>(true, {ilog2(2 * adr - 1), cmd.lz_ofs()}, {5, cmd.len & 0x1f}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 2 * 8 == ret.bit_length());
  write16(ret.out, 0, input.size());
  return ret.out;
}

} // namespace sfc_comp
