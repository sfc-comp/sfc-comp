#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> slap_stick_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0x10000);

  enum tag { uncomp, lz };

  static constexpr size_t pad = 0x11;

  std::vector<uint8_t> input(in.size() + pad);
  std::copy(in.begin(), in.end(), input.begin() + pad);
  for (size_t i = 0; i < pad; ++i) input[i] = 0x20;

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.c<0>(0x11);

  for (size_t i = input.size(); i-- > pad; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    dp.update(i, 2, 0x11, lz_helper.find(i, 0x100, 2), c0, 13, lz);
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(2);
  size_t adr = pad;
  for (const auto& cmd : dp.optimal_path(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<bnh>({9, size_t(0x100 | input[adr])}); break;
    case lz: ret.write<bnh, bnh>({9, (cmd.lz_ofs() - pad - 0x11) & 0xff}, {4, cmd.len - 2}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, in.size());

  assert(adr - pad == in.size());
  assert(dp.optimal_cost(pad) + 2 * 8 == ret.bit_length());

  return ret.out;
}

std::vector<uint8_t> sotsugyou_bangai_hen_comp(std::span<const uint8_t> input) {
  auto ret = slap_stick_comp(input);
  for (size_t i = 0; i < ret.size() - 2; ++i) ret[i] = ret[i + 2];
  ret.resize(ret.size() - 2);
  return ret;
}

} // namespace sfc_comp
