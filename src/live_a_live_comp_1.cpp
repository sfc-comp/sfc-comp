#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> live_a_live_comp_1(std::span<const uint8_t> input) {
  enum tag {
    uncomp,
    rle8z, rle8,
    rle16z, rle16,
    rle24,
    common_lo16, common_lo24, common_lo32,
    inc8, dec8, add8, add16,
    lzl, lzs, lz8
  };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x113); auto c1 = dp.c<1>(0xf0);
  auto c0_2 = dp.c<0, 2>(0x204); auto c1_2 = dp.c<1, 2>(0x206);
  auto c0_3 = dp.c<0, 3>(0x303); auto c1_3 = dp.c<1, 3>(0x306);
  auto c1_4 = dp.c<1, 4>(0x404);

  size_t rlen8 = 0, rlen16 = 0, rlen24 = 0;
  size_t rlen8d = 0; std::array<size_t, 2> rlen16d = {};
  std::array<size_t, 2> c16 = {};
  std::array<size_t, 3> c24 = {};
  std::array<size_t, 4> c32 = {};
  std::array<size_t, 16> lz8s = {};

  if (input.size() > 0) lz_helper.reset(input.size() - 1);
  for (size_t i = input.size(); i-- > 0; ) {
    if (i > 0) lz_helper.reset(i - 1);
    dp.update(i, 1, 0xf0, c1, 1, uncomp);

    rlen8 = encode::run_length_r(input, i, rlen8);
    if (!(input[i] & 0xf0)) {
      dp.update(i, 3, 0x12, rlen8, c0, 2, rle8z);
      dp.update(i, 0x13, 0x103, rlen8, c0, 3, rle8);
    } else {
      dp.update(i, 4, 0x103, rlen8, c0, 3, rle8);
    }

    rlen16 = encode::run_length16_r(input, i, rlen16);
    if (rlen16 >= 2) {
      if (((input[i] | input[i + 1]) & 0xf0) == 0) {
        dp.update(i, 4, 0x202, rlen16, c0_2, 3, rle16z);
      } else {
        dp.update(i, 4, 0x202, rlen16, c0_2, 4, rle16);
      }
    }

    rlen24 = encode::run_length24_r(input, i, rlen24);
    dp.update(i, 6, 0x303, rlen24, c0_3, 5, rle24);

    c16[i % 2] = encode::common_lo16_r(input, i, c16[i % 2]);
    dp.update(i, 8, 0x206, c16[i % 2], c1_2, 3, common_lo16);

    c24[i % 3] = encode::common_lo24_16_r(input, i, c24[i % 3]);
    dp.update(i, 9, 0x306, c24[i % 3], c1_3, 4, common_lo24);

    c32[i % 4] = encode::common_lo32_24_r(input, i, c32[i % 4]);
    dp.update(i, 8, 0x404, c32[i % 4], c1_4, 5, common_lo32);

    rlen8d = encode::run_length_delta_r(input, i, rlen8d);
    if (rlen8d >= 2) {
      const uint8_t delta = input[i + 1] - input[i];
      if (delta == 1 || delta == 0xff) {
        dp.update(i, 4, 0x103, rlen8d, c0, 3, (delta == 1) ? inc8 : dec8);
      } else {
        dp.update(i, 5, 0x104, rlen8d, c0, 4, add8);
      }
    }
    rlen16d[i % 2] = encode::run_length16_delta_r(input, i, rlen16d[i % 2]);
    if (rlen16d[i % 2] >= 4) {
      const uint16_t delta = read16(input, i + 2) - read16(input, i);
      if (((delta + 0x80) & 0xffff) < 0x100) {
        dp.update(i, 6, 0x204, rlen16d[i % 2], c0_2, 5, add16);
      }
    }
    dp.update(i, 4, 0x13, lz_helper.find(i, 0x1000, 4), c0, 3, lzl);
    dp.update(i, 0x14, 0x113, lz_helper.find(i, 0x100, 0x14), c0, 3, lzs);

    for (size_t k = 0; k < 0x10; ++k) lz8s[k] = encode::lz_dist_r(input, i, 8 * (k + 1), lz8s[k]);
    const size_t best_k = std::max_element(lz8s.begin(), lz8s.end()) - lz8s.begin();
    dp.update(i, 3, 0x12, {i - 8 * (best_k + 1), lz8s[best_k]}, c0, 2, lz8);

    c0.update(i); c1.update(i); c0_2.update(i); c1_2.update(i);
    c0_3.update(i); c1_3.update(i); c1_4.update(i);
  }

  using namespace data_type;
  writer ret; ret.write<d8>(1);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
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
      ret.write<d8, d8, d16, d8>(0xfb, (cmd.len - 6) >> 1,
        read16(input, adr), input[adr + 2] - input[adr]
      );
    } break;
    case lzl: {
      const size_t d = adr - cmd.lz_ofs(); assert(d >= 2);
      ret.write<d8, d16>(0xfc, (d - 1) | (cmd.len - 4) << 12);
    } break;
    case lzs: {
      const size_t d = adr - cmd.lz_ofs(); assert(d >= 2);
      ret.write<d8, d8, d8>(0xfd, d - 1, cmd.len - 20);
    } break;
    case lz8: {
      const size_t d = adr - cmd.lz_ofs();
      assert(!(d & 7));
      ret.write<d8, d8>(0xfe, (d - 8) << 1 | (cmd.len - 3));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<d8>(0xff);
  assert(dp.optimal_cost() + 2 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
