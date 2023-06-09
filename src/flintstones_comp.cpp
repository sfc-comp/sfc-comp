#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> flintstones_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  static constexpr auto dists_s = std::to_array<size_t>({
    0x0001, 0x0002, 0x0003, 0x0004, 0x0006, 0x0008, 0x000c, 0x0010,
    0x0018, 0x0020, 0x0030, 0x0040, 0x0060, 0x0080, 0x00c0, 0x0100
  });

  static constexpr auto dists_l = std::to_array<size_t>({
    0x0001, 0x0002, 0x0003, 0x0004, 0x0006, 0x0008, 0x000c, 0x0010,
    0x0018, 0x0020, 0x0030, 0x0040, 0x0060, 0x0080, 0x00a0, 0x00c0,
    0x00e0, 0x0100, 0x0140, 0x0180, 0x0200, 0x0300, 0x0400, 0x0600,
    0x0800, 0x0c00, 0x1000, 0x1800, 0x2000, 0x2800, 0x3000, 0x4000
  });

  enum tag {
    uncomp,
    lzs, lzm,
    lzds, lzdl,
    lzl
  };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(9 + 0x3f);

  std::array<size_t, dists_s.size()> res_lz_ds = {};
  std::array<size_t, dists_l.size()> res_lz_dl = {};

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 9, uncomp);

    dp.update(i, 2, 4, lz_helper.find(i, 0x40, 2), c0, 9, lzs);
    dp.update(i, 3, 8, lz_helper.find(i, 0x1000, 3), c0, 18, lzm);

    for (size_t k = 0; k < dists_s.size(); ++k) res_lz_ds[k] = encode::lz_dist_r(input, i, dists_s[k], res_lz_ds[k]);
    const size_t ks = std::max_element(res_lz_ds.begin(), res_lz_ds.end()) - res_lz_ds.begin();
    dp.update(i, 5, 12, {ks, res_lz_ds[ks]}, c0, 13, lzds);

    for (size_t k = 0; k < dists_l.size(); ++k) res_lz_dl[k] = encode::lz_dist_r(input, i, dists_l[k], res_lz_dl[k]);
    const size_t kl = std::max_element(res_lz_dl.begin(), res_lz_dl.end()) - res_lz_dl.begin();
    dp.update(i, 9, 9 + 0x3f, {kl, res_lz_dl[kl]}, c0, 18, lzdl);

    auto res_lzl = lz_helper.find_closest(i, 0x4000, 4, 0x43);
    if ((i - res_lzl.ofs) > 0x3e00 && res_lzl.len == 0x43) res_lzl.len -= 1;
    dp.update(i, 4, 4 + 0x3f, res_lzl, c0, 27, lzl);

    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(3);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: {
      ret.write<b1, d8>(false, input[adr]);
    } break;
    case lzs: {
      ret.write<b1, d8>(true, 0x00 | (cmd.len - 2) << 6 | (d - 1));
    } break;
    case lzm: {
      ret.write<b1, d8>(true, 0xc0 | (cmd.len - 3) << 3 | (d - 1) >> 9);
      ret.write<b1, d8>(((d - 1) >> 8) & 1, (d - 1) & 0x00ff);
    } break;
    case lzds: {
      const size_t k = cmd.arg;
      ret.write<b1, d8>(true, 0xf0 | (cmd.len - 5));
      ret.write<bnh>({4, k});
    } break;
    case lzdl: {
      const size_t k = cmd.arg;
      ret.write<b1, d8>(true, 0xf8 | k >> 3);
      ret.write<b1, d8>((k >> 2) & 1, (k & 0x03) << 6 | (cmd.len - 9));
    } break;
    case lzl: {
      assert(!(cmd.len == 0x43 && (d - 1) >= 0x3e00));
      ret.write<b1, d8>(true, 0xfc | (cmd.len - 4) >> 4);
      ret.write<b1, d8>(((cmd.len - 4) >> 3) & 1, ((cmd.len - 4) & 0x07) << 5 | (d - 1) >> 9);
      ret.write<b1, d8>(((d - 1) >> 8) & 1, (d - 1) & 0x00ff);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1, d8, b1, d8>(true, 0xff, true, 0xff);
  assert(adr == input.size());
  assert(dp.optimal_cost() + 18 + 3 * 8 == ret.bit_length());

  ret[0] = 0x9a;
  write16(ret.out, 1, input.size());
  return ret.out;
}

} // namespace sfc_comp
