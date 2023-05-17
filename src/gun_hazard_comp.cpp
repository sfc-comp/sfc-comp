#include <queue>

#include "image.hpp"
#include "lzss.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> gun_hazard_comp_1(std::span<const uint8_t> input, const uint8_t header_val) {
  check_divisibility(input.size(), 32);
  check_size(input.size(), 0x0020, 0x10000);

  using namespace data_type;
  writer_b4_h ret(3);
  const auto indexed_256 = snes4bpp::to_indexed256_8_1(input);
  for (size_t i = 0; i < indexed_256.size(); ) {
    const auto v = indexed_256[i];
    ret.write<b4>(v);
    size_t len = encode::run_length(indexed_256, i, 0);
    i += len;
    for (--len; len > 0; ) {
      size_t l = std::min<size_t>(0x10, len);
      ret.write<b4, b4>(v, l - 1),
      len -= l;
    }
  }
  ret[0] = header_val | 0x01;
  write16(ret.out, 1, input.size());
  return ret.out;
}

std::vector<uint8_t> gun_hazard_comp_2(std::span<const uint8_t> input, const uint8_t header_val) {
  check_divisibility(input.size(), 0x20);
  check_size(input.size(), 0x0020, 0x10000);

  std::array<std::vector<size_t>, 3> zeros;
  for (size_t zi = 0; zi < zeros.size(); ++zi) {
    zeros[zi] = std::vector<size_t>(1 << 16);
  }
  std::array<size_t, 16> freq = {};

  auto indexed_256 = snes4bpp::to_indexed256_8_1(input);
  for (size_t i = 0; i < indexed_256.size(); i += 32) {
    std::array<uint16_t, 4> pal_flags = {};
    for (size_t j = 0; j < 4; ++j) {
      uint16_t pal_flag = 0;
      for (size_t k = 0; k < 8; ++k) {
        const auto v = indexed_256[i + j * 8 + k];
        freq[v] += 1;
        pal_flag |= 1 << v;
      }
      pal_flags[j] = pal_flag, zeros[0][pal_flag] += 1;
      if (j & 1) {
        pal_flags[j] |= pal_flags[j - 1], zeros[1][pal_flags[j]] += 1;
        if (j & 2) {
          pal_flags[j] |= pal_flags[j - 2], zeros[2][pal_flags[j]] += 1;
        }
      }
    }
  }

  const auto zeta = [](std::vector<size_t>& cumu) {
    for (size_t lg = 0; lg < 16; ++lg) {
      const size_t m = 2 << lg, mh = m >> 1;
      for (size_t i = 0; i < 0x10000; i += m) {
        for (size_t j = 0; j < mh; ++j) cumu[i + j + mh] += cumu[i + j];
      }
    }
  };

  for (size_t zi = 0; zi < zeros.size(); ++zi) zeta(zeros[zi]);

  using perm_type = std::array<size_t, 16>;

  const perm_type iden_perm = []{
    perm_type ret;
    std::iota(ret.begin(), ret.end(), 0);
    return ret;
  }();

  perm_type best_perm = iden_perm;
  size_t best_cost = std::numeric_limits<size_t>::max();
  size_t best_conf = 0;

  const auto update = [&](const perm_type& perm) -> bool {
    size_t bitplane_flags[4] = {};
    for (size_t i = 0; i < perm.size(); ++i) {
      for (size_t j = 0; j < 4; ++j) {
        if (!((i >> j) & 1)) bitplane_flags[j] |= 1 << perm[i];
      }
    }
    const bool iden = (perm == iden_perm);
    static constexpr auto lens = std::to_array<size_t>({1, 2, 4});
    static constexpr auto numers = std::to_array<size_t>({36, 34, 33});
    bool updated = false;
    for (size_t zi = 0; zi < 3; ++zi) {
      size_t z = 0;
      for (auto f : bitplane_flags) z += zeros[zi][f];
      size_t cost = input.size() * numers[zi] / 32 - z * lens[zi] + (iden ? 0 : 8);
      if (cost < best_cost) {
        best_cost = cost;
        best_perm = perm;
        best_conf = lens[zi] | (iden ? 0x00 : 0x80);
        updated = true;
      }
    }
    return updated;
  };

  const auto initial_perm = [&]{
    perm_type perm = iden_perm;
    std::sort(perm.begin(), perm.end(), [&](const size_t i, const size_t j) {
      return freq[i] > freq[j];
    });
    // [TODO] might not be useful
    static constexpr auto pos = perm_type({
      0, 1, 2, 4, 8, 3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15
    });
    perm_type ret;
    for (size_t i = 0; i < ret.size(); ++i) ret[pos[i]] = perm[i];
    return ret;
  }();

  // [TODO] Find a better method.
  update(iden_perm);
  update(initial_perm);

  while (true) {
    perm_type perm = best_perm;
    bool updated = false;
    for (size_t i = 0; i < perm.size(); ++i) {
      for (size_t j = i + 1; j < perm.size(); ++j) {
        std::swap(perm[i], perm[j]);
        for (size_t k = i; k < perm.size(); ++k) {
          for (size_t l = k; l < perm.size(); ++l) {
            std::swap(perm[k], perm[l]);
            if (update(perm)) updated = true;
            std::swap(perm[k], perm[l]);
          }
        }
        std::swap(perm[i], perm[j]);
      }
    }
    if (!updated) break;
  }

  using namespace data_type;
  writer_b8_l ret(4);
  if (best_conf & 0x80) {
    for (size_t i = 0; i < 16; i += 2) ret.write<d8>(best_perm[i] << 4 | best_perm[i + 1]);
    perm_type iperm;
    for (size_t i = 0; i < best_perm.size(); ++i) iperm[best_perm[i]] = i;
    for (size_t i = 0; i < indexed_256.size(); ++i) indexed_256[i] = iperm[indexed_256[i]];
  }

  const size_t block_size = best_conf & 0x7f;
  const auto permuted_4bpp = snes4bpp::from_indexed256_8_1(indexed_256);
  for (size_t i = 0; i < permuted_4bpp.size(); i += 0x20) {
    for (const size_t offset : {0x00, 0x10, 0x01, 0x11}) {
      for (size_t o = offset; o < offset + 0x10; o += block_size * 2) {
        bool zero = true;
        for (size_t j = 0; j < block_size; ++j) {
          if (permuted_4bpp[i + o + 2 * j] != 0) {
            zero = false;
            break;
          }
        }
        ret.write<b1>(!zero);
        if (zero) continue;
        for (size_t j = 0; j < block_size; ++j) {
          ret.write<d8>(permuted_4bpp[i + o + 2 * j]);
        }
      }
    }
  }
  assert(ret.size() == best_cost + 4);

  ret[0] = header_val | 0x02;
  write16(ret.out, 1, input.size());
  ret[3] = best_conf;
  return ret.out;
}

std::vector<uint8_t> gun_hazard_comp_3(std::span<const uint8_t> input, const uint8_t header_val) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss<writer_b8_l>(
    input, 0, [](std::span<const uint8_t>) {},
    0x0fff, 3, 0x12,
    3, true,
    [&](size_t i, size_t o, size_t l) {
      size_t d = (i - o);
      return (d & 0x0f00) << 4 | (l - 3) << 8 | (d & 0x00ff);
    }
  );
  ret[0] = header_val | 0x03;
  write16(ret, 1, input.size());
  return ret;
}

std::vector<uint8_t> gun_hazard_comp_4(std::span<const uint8_t> input, const uint8_t header_val) {
  check_size(input.size(), 1, 0x10000);

  enum tag { uncomp, lzs, lzl, lzll };

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(256);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    dp.update(i, 3, 6, lz_helper.find(i, 0xff, 3), c0, 12, lzs);
    const auto res_lzl = lz_helper.find(i, 0x1fff, 3);
    dp.update(i, 3, 9, res_lzl, c0, 18, lzl);
    dp.update(i, 1, 256, res_lzl, c0, 26, lzll);
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_l ret(3);

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: ret.write<b1, d8>(true, input[adr]); break;
    case lzs: ret.write<bnh, d8>({4, cmd.len - 3}, d); break;
    case lzl: ret.write<b1, b1, d8, d8>(false, true, d & 0xff,  (d >> 8) << 3 | (cmd.len - 2)); break;
    case lzll: ret.write<b1, b1, d8, d8, d8>(false, true, d & 0xff, (d >> 8) << 3 | 0, cmd.len - 1); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 3 * 8 == ret.bit_length());

  ret[0] = header_val | 0x04;
  write16(ret.out, 1, input.size());
  return ret.out;
}

std::vector<uint8_t> gun_hazard_comp_5(std::span<const uint8_t> input, const uint8_t header_val) {
  check_divisibility(input.size(), 0x20);
  auto ret = gun_hazard_comp_3(snes4bpp::to_indexed16_h_8_1(input), header_val);
  ret[0] = header_val | 0x05;
  return ret;
}

std::vector<uint8_t> gun_hazard_comp_6(std::span<const uint8_t> input, const uint8_t header_val) {
  check_divisibility(input.size(), 0x20);
  auto ret = gun_hazard_comp_4(snes4bpp::to_indexed16_h_8_1(input), header_val);
  ret[0] = header_val | 0x06;
  return ret;
}

std::vector<uint8_t> gun_hazard_comp_core(
    std::span<const uint8_t> input, const uint8_t header_val, const size_t version) {
  check_size(input.size(), 1, 0x10000);

  const bool maybe_4bpp = (input.size() % 32 == 0);

  std::vector<uint8_t> best(input.size() + 3);
  std::copy(input.begin(), input.end(), best.begin() + 3);
  best[0] = header_val | 0x00;
  write16(best, 1, input.size());

  if (maybe_4bpp) {
    if (auto res = gun_hazard_comp_1(input, header_val); res.size() < best.size()) {
      best = std::move(res);
    }
    if (auto res = gun_hazard_comp_2(input, header_val); res.size() < best.size()) {
      best = std::move(res);
    }
  }

  if (auto res = gun_hazard_comp_3(input, header_val); res.size() < best.size()) {
    best = std::move(res);
  }
  if (version >= 1) {
    if (auto res = gun_hazard_comp_4(input, header_val); res.size() < best.size()) {
      best = std::move(res);
    }
    if (version >= 2 && maybe_4bpp) {
      if (auto res = gun_hazard_comp_5(input, header_val); res.size() < best.size()) {
        best = std::move(res);
      }
      if (auto res = gun_hazard_comp_6(input, header_val); res.size() < best.size()) {
        best = std::move(res);
      }
    }
  }
  return best;
}

} // namespace

std::vector<uint8_t> assault_suits_valken_comp(std::span<const uint8_t> input) {
  // Note: The value 0xa0 can be 0x00 for some cases.
  return gun_hazard_comp_core(input, 0xa0, 0);
}

std::vector<uint8_t> brandish_comp(std::span<const uint8_t> input) {
  return gun_hazard_comp_core(input, 0xa0, 1);
}

std::vector<uint8_t> gun_hazard_comp(std::span<const uint8_t> input) {
  return gun_hazard_comp_core(input, 0x20, 2);
}

} // namespace sfc_comp
