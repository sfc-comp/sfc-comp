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

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x3ff);
  auto c1 = dp.c<1>(0x20);
  auto c1_2 = dp.c<1, 2>(0x22);

  size_t rlen = 0, lo8_len = 0, hi8_len = 0;
  std::array<size_t, 2> lo16_lens = {};
  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 0x20, c1, 1, uncomp);
    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, 1, 0x20, rlen, c0, 1, rle0);
      dp.update(i, 0x21, 0x3ff, rlen, c0, 2, rle0l);
    } else if (input[i] == 0xff) {
      dp.update(i, 1, 8, rlen, c0, 1, rle1);
      dp.update(i, 9, 0x11, rlen, c0, 2, rle);
      dp.update(i, 0x12, 0x3ff, rlen, c0, 3, rlel);
    } else {
      dp.update(i, 2, 0x11, rlen, c0, 2, rle);
      dp.update(i, 0x12, 0x3ff, rlen, c0, 3, rlel);
    }

    lo16_lens[i & 1] = encode::common_lo16_r(input, i, lo16_lens[i & 1]);
    if (input[i] == 0) dp.update(i, 2, 0x20, lo16_lens[i & 1], c1_2, 1, common_lo16_0);
    dp.update(i, 4, 0x22, lo16_lens[i & 1], c1_2, 2, common_lo16);

    hi8_len = encode::common_hi8_r(input, i, hi8_len);
    if ((input[i] & 0xf0) == 0) {
      dp.update(i, 2, 0x20, hi8_len, c1_2, 1, common_hi8_0);
    } else if ((input[i] & 0xf0) == 0xf0) {
      dp.update(i, 2, 0x20, hi8_len, c1_2, 1, common_hi8_f);
    }
    lo8_len = encode::common_lo8_r(input, i, lo8_len);
    if ((input[i] & 0x0f) == 0) {
      dp.update(i, 2, 0x20, lo8_len, c1_2, 1, common_lo8_0);
    } else if ((input[i] & 0x0f) == 0x0f) {
      dp.update(i, 2, 0x20, lo8_len, c1_2, 1, common_lo8_f);
    }
    dp.update(i, 2, 0x11, lz_helper.find(i, 0x400, 2), c0, 2, lz);

    c0.update(i); c1.update(i); c1_2.update(i);
  }

  using namespace data_type;
  writer ret(2);
  size_t adr = pad;
  for (const auto& cmd : dp.optimal_path(pad)) {
    switch (cmd.type) {
    case lz: {
      ret.write<d16b>((cmd.len - 2) << 10 | ((adr - cmd.lz_ofs()) & 0x03ff));
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
  assert(dp.optimal_cost(pad) + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
