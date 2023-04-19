#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> live_a_live_comp_1(std::span<const uint8_t> input) {
  enum CompType {
    uncomp,
    rle8z, rle8,
    rle16z, rle16,
    rle24,
    common_lo16, common_lo24, common_lo32,
    inc8, dec8, add8, add16,
    lzl, lzs, lz8
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  size_t rlen8 = 0, rlen16 = 0, rlen24 = 0;
  encode::rle_data c16[2] = {{}, {}}, c24[3] = {{}, {}, {}}, c32[4] = {{}, {}, {}, {}};
  encode::rle_data rlen8i = {}, rlen16i[2] = {{}, {}};
  std::array<size_t, 16> lz8s = {};
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 0xf0, Linear<1, 1>(), uncomp);

    rlen8 = encode::run_length(input, i, rlen8);
    if (!(input[i] & 0xf0)) {
      dp.update(i, 3, 0x12, rlen8, Constant<2>(), rle8z);
      dp.update(i, 0x13, 0x103, rlen8, Constant<3>(), rle8);
    } else {
      dp.update(i, 4, 0x103, rlen8, Constant<3>(), rle8);
    }

    if (i + 1 < input.size()) {
      rlen16 = encode::run_length16(input, i, rlen16);
      if (!((input[i] | input[i + 1]) & 0xf0)) {
        dp.update_k<2>(i, 4, 0x202, rlen16, Constant<3>(), rle16z);
      } else {
        dp.update_k<2>(i, 4, 0x202, rlen16, Constant<4>(), rle16);
      }
    }

    rlen24 = encode::run_length24(input, i, rlen24);
    dp.update_k<3>(i, 6, 0x0303, rlen24, Constant<5>(), rle24);

    size_t rem2 = i % 2;
    c16[rem2] = encode::common_lo16_hint(input, i, c16[rem2].len);
    dp.update_k<2>(i, 8, 0x206, c16[rem2].len, LinearQ<1, 6 + 1, 2>(), common_lo16);

    size_t rem3 = i % 3;
    c24[rem3] = encode::common_lo24_16_hint(input, i, c24[rem3].len);
    dp.update_k<3>(i, 9, 0x306, c24[rem3].len, LinearQ<1, 12 + 2, 3>(), common_lo24);

    size_t rem4 = i % 4;
    c32[rem4] = encode::common_lo32_24_hint(input, i, c32[rem4].len);
    dp.update_k<4>(i, 8, 0x404, c32[rem4].len, LinearQ<1, 20 + 3, 4>(), common_lo32);

    if (i + 1 < input.size()) {
      rlen8i = encode::run_length_delta(input, i, rlen8i);
      size_t delta = rlen8i.v;
      if (delta == 1 || delta == 0xff) {
        dp.update(i, 4, 0x103, rlen8i.len, Constant<3>(), (delta == 1) ? inc8 : dec8);
      } else {
        dp.update(i, 5, 0x104, rlen8i.len, Constant<4>(), add8);
      }
    }

    if (i + 3 < input.size()) {
      rlen16i[rem2] = encode::run_length16_delta(input, i, rlen16i[rem2]);
      size_t delta = rlen16i[rem2].v;
      if (((delta + 0x80) & 0xffff) < 0x100) {
        dp.update_k<2>(i, 6, 0x204, rlen16i[rem2].len, Constant<5>(), add16);
      }
    }

    auto res_lzl = lz_helper.find_best(i, 0x1000);
    dp.update_lz(i, 4, 0x13, res_lzl, Constant<3>(), lzl);

    auto res_lzs = lz_helper.find_best(i, 0x100);
    dp.update_lz(i, 0x14, 0x113, res_lzs, Constant<3>(), lzs);

    for (size_t k = 0; k < 0x10; ++k) lz8s[k] = encode::lz_dist(input, i, 8 * (k + 1), lz8s[k]);
    const size_t best_k = std::max_element(lz8s.begin(), lz8s.end()) - lz8s.begin();
    dp.update_lz(i, 3, 0x12, {ptrdiff_t(i - 8 * (best_k + 1)), lz8s[best_k]}, Constant<2>(), lz8);

    // offset should be >= 2.
    if (i >= 1) lz_helper.add_element(i - 1);
  }

  using namespace data_type;
  writer ret; ret.write<d8>(1);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      ret.write<d8, d8n>(cmd.len - 1, {cmd.len, &input[adr]});
    } break;
    case rle8z: {
      ret.write<d8, d8>(0xf0, input[adr] << 4 | (cmd.len - 3));
    } break;
    case rle8: {
      ret.write<d8, d8, d8>(0xf1, cmd.len - 4, input[adr]);
    } break;
    case rle16z: {
      ret.write<d8, d8, d8>(0xf2, (cmd.len - 4) >> 1, input[adr + 1] << 4 | input[adr]);
    } break;
    case rle16: {
      ret.write<d8, d8, d8, d8>(0xf3, (cmd.len - 4) >> 1, input[adr], input[adr + 1]);
    } break;
    case rle24: {
      ret.write<d8, d8, d8, d8, d8>(0xf4, (cmd.len - 6) / 3, input[adr], input[adr + 1], input[adr + 2]);
    } break;
    case common_lo16: {
      ret.write<d8, d8, d8, d8nk>(0xf5, (cmd.len - 8) >> 1, input[adr], {cmd.len, &input[adr + 1], 2});
    } break;
    case common_lo24: {
      ret.write<d8, d8, d8, d8, d8nk>(
        0xf6, (cmd.len - 9) / 3, input[adr], input[adr + 1], {cmd.len, &input[adr + 2], 3}
      );
    } break;
    case common_lo32: {
      ret.write<d8, d8, d8, d8, d8, d8nk>(
        0xf7, (cmd.len - 8) / 4, input[adr], input[adr + 1], input[adr + 2],
        {cmd.len, &input[adr + 3], 4}
      );
    } break;
    case inc8: {
      ret.write<d8, d8, d8>(0xf8, cmd.len - 4, input[adr]);
    } break;
    case dec8: {
      ret.write<d8, d8, d8>(0xf9, cmd.len - 4, input[adr]);
    } break;
    case add8: {
      ret.write<d8, d8, d8, d8>(0xfa, cmd.len - 5, input[adr], input[adr + 1] - input[adr]);
    } break;
    case add16: {
      ret.write<d8, d8, d16, d8>(
        0xfb, (cmd.len - 6) >> 1,
        read16(input, adr), input[adr + 2] - input[adr]
      );
    } break;
    case lzl: {
      size_t d = adr - cmd.lz_ofs; assert(d >= 2);
      ret.write<d8, d16>(0xfc, (d - 1) | (cmd.len - 4) << 12);
    } break;
    case lzs: {
      size_t d = adr - cmd.lz_ofs; assert(d >= 2);
      ret.write<d8, d8, d8>(0xfd, d - 1, cmd.len - 20);
    } break;
    case lz8: {
      size_t d = adr - cmd.lz_ofs;
      assert(!(d & 7));
      ret.write<d8, d8>(0xfe, (d - 8) << 1 | (cmd.len - 3));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d8>(0xff);
  assert(dp.total_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
