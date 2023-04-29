#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> konami_comp_core(std::span<const uint8_t> in, const bool use_old_version) {
  check_size(in.size(), 0, 0x10000);
  enum CompType {
    uncomp, rle0, rle0l, rle, lz, common16
  };
  static constexpr size_t pad = 0x21;
  std::vector<uint8_t> input(in.size() + pad);
  std::copy(in.begin(), in.end(), input.begin() + pad);

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);
  for (size_t i = 0; i < pad; ++i) dp[i + 1].cost = 0;

  size_t rlen = 0;
  for (size_t i = pad; i < input.size(); ++i) {
    dp.update(i, 1, 31, Linear<1, 1>(), uncomp);
    rlen = encode::run_length(input, i, rlen);
    if (input[i] != 0) {
      dp.update(i, 2, 0x21, rlen, Constant<2>(), rle);
    } else {
      if (use_old_version) {
        dp.update(i, 2, 0x21, rlen, Constant<1>(), rle0);
      } else {
        dp.update(i, 2, 0x20, rlen, Constant<1>(), rle0);
        dp.update(i, 0x21, 0x101, rlen, Constant<2>(), rle0l);
      }
    }
    auto res_lz = lz_helper.find_best(i, 0x400);
    dp.update_lz(i, 2, 0x21, res_lz, Constant<2>(), lz);
    if (input[i] == 0) {
      auto common16_len = encode::common_lo16(input, i, 0x42).len;
      dp.update_k<2>(i, 4, 0x42, common16_len, LinearQ<1, 2, 2>(), common16);
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = pad;
  for (const auto cmd : dp.commands(pad)) {
    size_t d = (cmd.lz_ofs - pad - 0x21) & 0x03ff;
    switch (cmd.type) {
    case lz: ret.write<d16b>((cmd.len - 2) << 10 | d); break;
    case uncomp: ret.write<d8, d8n>(0x80 + cmd.len, {cmd.len, &input[adr]}); break;
    case common16: ret.write<d8, d8nk>(
      0xa0 + ((cmd.len - 4) >> 1), {cmd.len, &input[adr + 1], 2}
    ); break;
    case rle: ret.write<d8, d8>(0xc0 + cmd.len - 2, input[adr]); break;
    case rle0: ret.write<d8>(0xe0 + cmd.len - 2); break;
    case rle0l: ret.write<d8, d8>(0xff, cmd.len - 2); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, ret.size());
  assert(dp.total_cost() + 2 == ret.size());
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
