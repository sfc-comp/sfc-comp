#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> konami_comp_core(std::span<const uint8_t> in, const bool use_old_version) {
  check_size(in.size(), 0, 0x10000);
  enum tag { uncomp, rle0, rle0l, rle, lz, common16 };
  static constexpr size_t pad = 0x21;
  std::vector<uint8_t> input(in.size() + pad);
  std::ranges::copy(in, input.begin() + pad);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x101);
  auto c1 = dp.c<1>(31);
  auto c1_2 = dp.c<1, 2>(0x42);

  size_t rlen = 0; std::array<size_t, 2> lo16_lens = {};
  for (size_t i = input.size(); i-- > pad; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 31, c1, 1, uncomp);
    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] != 0) {
      dp.update(i, 2, 0x21, rlen, c0, 2, rle);
    } else {
      if (use_old_version) {
        dp.update(i, 2, 0x21, rlen, c0, 1, rle0);
      } else {
        dp.update(i, 2, 0x20, rlen, c0, 1, rle0);
        dp.update(i, 0x21, 0x101, rlen, c0, 2, rle0l);
      }
    }
    dp.update(i, 2, 0x21, lz_helper.find(i, 0x400, 2), c0, 2, lz);
    lo16_lens[i & 1] = encode::common_lo16_r(input, i, lo16_lens[i & 1]);
    if (input[i] == 0) dp.update(i, 4, 0x42, lo16_lens[i & 1], c1_2, 1, common16);
    c0.update(i); c1.update(i); c1_2.update(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = pad;
  for (const auto& cmd : dp.optimal_path(pad)) {
    const size_t d = (cmd.lz_ofs() - pad - 0x21) & 0x03ff;
    switch (cmd.type) {
    case lz: ret.write<d16b>((cmd.len - 2) << 10 | d); break;
    case uncomp: ret.write<d8, d8n>(0x80 + cmd.len, {cmd.len, &input[adr]}); break;
    case common16: ret.write<d8, d8nk>(0xa0 + ((cmd.len - 4) >> 1),
                                       {cmd.len, &input[adr + 1], 2}); break;
    case rle: ret.write<d8, d8>(0xc0 + cmd.len - 2, input[adr]); break;
    case rle0: ret.write<d8>(0xe0 + cmd.len - 2); break;
    case rle0l: ret.write<d8, d8>(0xff, cmd.len - 2); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, ret.size());
  assert(dp.optimal_cost(pad) + 2 == ret.size());
  assert(adr - pad == in.size());
  return ret.out;
}

} // namespace

std::vector<uint8_t> konami_comp_1(std::span<const uint8_t> input) {
  return konami_comp_core(input, true);
}

std::vector<uint8_t> konami_comp_2(std::span<const uint8_t> input) {
  return konami_comp_core(input, false);
}

} // namespace sfc_comp
