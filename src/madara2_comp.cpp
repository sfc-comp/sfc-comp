#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> madara2_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0x10000);
  enum tag {
    uncomp,
    rle0, rle0l,
    rle1,
    rle, rlel,
    lz,
    common_lo8_0, common_lo8_f,
    common_hi8_0, common_hi8_f,
    common_lo16_0, common_lo16
  };

  static constexpr size_t pad = 0x11;
  std::vector<uint8_t> input(in.size() + pad);
  std::ranges::copy(in, input.begin() + pad);

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size(), pad);
  for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);

  size_t rlen = 0;
  for (size_t i = pad; i < input.size(); ++i) {
    dp.update(i, 1, 0x20, Linear<1, 1>(), uncomp);
    rlen = encode::run_length(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, 1, 0x20, rlen, Constant<1>(), rle0);
      dp.update(i, 0x21, 0x3ff, rlen, Constant<2>(), rle0l);
    } else if (input[i] == 0xff) {
      dp.update(i, 1, 8, rlen, Constant<1>(), rle1);
      dp.update(i, 9, 0x11, rlen, Constant<2>(), rle);
      dp.update(i, 0x12, 0x3ff, rlen, Constant<3>(), rlel);
    } else {
      dp.update(i, 2, 0x11, rlen, Constant<2>(), rle);
      dp.update(i, 0x12, 0x3ff, rlen, Constant<3>(), rlel);
    }

    auto common_lo16_len = encode::common_lo16(input, i, 0x22).len;
    if (input[i] == 0) {
      dp.update_k<2>(i, 2, 0x20, common_lo16_len, LinearQ<1, 2, 2>(), common_lo16_0);
    }
    dp.update_k<2>(i, 4, 0x22, common_lo16_len, LinearQ<1, 4, 2>(), common_lo16);

    auto common_hi8_res = encode::common_hi8(input, i, 0x20);
    if (common_hi8_res.v == 0) {
      dp.update_k<2>(i, 2, 0x20, common_hi8_res.len,
        LinearQ<1, 2, 2>(), common_hi8_0);
    } else if (common_hi8_res.v == 0xf0) {
      dp.update_k<2>(i, 2, 0x20, common_hi8_res.len,
        LinearQ<1, 2, 2>(), common_hi8_f);
    }
    auto common_lo8_res = encode::common_lo8(input, i, 0x20);
    if (common_lo8_res.v == 0) {
      dp.update_k<2>(i, 2, 0x20, common_lo8_res.len, LinearQ<1, 2, 2>(), common_lo8_0);
    } else if (common_lo8_res.v == 0x0f) {
      dp.update_k<2>(i, 2, 0x20, common_lo8_res.len, LinearQ<1, 2, 2>(), common_lo8_f);
    }
    auto res_lz = lz_helper.find(i, 0x400, 2);
    dp.update_lz(i, 2, 0x11, res_lz, Constant<2>(), lz);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = pad;
  for (const auto cmd : dp.commands(pad)) {
    switch (cmd.type) {
    case lz: {
      ret.write<d16b>((cmd.len - 2) << 10 | ((adr - cmd.lz_ofs) & 0x03ff));
    } break;
    case uncomp: {
      ret.write<d8, d8n>(0x40 + cmd.len - 1, {cmd.len, &input[adr]});
    } break;
    case rle0: {
      ret.write<d8>(0x60 + (cmd.len - 1));
    } break;
    case common_lo8_0: {
      ret.write<d8, d8n_hb>(0x80 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr]});
    } break;
    case common_hi8_0: {
      ret.write<d8, d8n_lb>(0x90 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr]});
    } break;
    case common_lo8_f: {
      ret.write<d8, d8n_hb>(0xa0 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr]});
    } break;
    case common_hi8_f: {
      ret.write<d8, d8n_lb>(0xb0 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr]});
    } break;
    case common_lo16_0: {
      ret.write<d8, d8nk>(0xc0 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr + 1], 2});
    } break;
    case rle: {
      ret.write<d8, d8>(0xd0 + (cmd.len - 2), input[adr]);
    } break;
    case common_lo16: {
      ret.write<d8, d8, d8nk>(0xe0 + ((cmd.len - 4) >> 1), input[adr], {cmd.len, &input[adr + 1], 2});
    } break;
    case rle1: {
      ret.write<d8>(0xf0 + (cmd.len - 1));
    } break;
    case rlel: {
      ret.write<d16b, d8>(0xf800 | cmd.len, input[adr]);
    } break;
    case rle0l: {
      ret.write<d16b>(0xfc00 | cmd.len);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, ret.size());
  assert(dp.optimal_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
