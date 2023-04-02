#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

namespace re {

constexpr size_t lzl_max_len = 0x22;
constexpr size_t lzl_max_ofs = 0x1fff;
constexpr size_t lzs_max_len = 10;
constexpr size_t lzs_max_ofs = 0x3f;
constexpr size_t max_len = 311;

constexpr size_t frees[max_len + 1] = {
  0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 21, 24, 24, 24,
  24, 24, 24, 24, 24, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24, 21, 21, 21, 27, 27,
  27, 27, 27, 24, 24, 24, 24, 24, 24, 24, 27, 27, 27, 24, 24, 24, 24, 24, 24, 30,
  30, 27, 27, 27, 27, 27, 27, 27, 27, 27, 24, 24, 24, 24, 27, 27, 27, 27, 27, 27,
  24, 24, 24, 30, 30, 30, 30, 30, 30, 30, 30, 27, 27, 27, 27, 27, 27, 27, 30, 30,
  30, 27, 27, 27, 27, 27, 27, 33, 33, 33, 33, 33, 30, 30, 30, 30, 30, 30, 30, 30,
  30, 27, 27, 27, 27, 30, 30, 30, 30, 30, 30, 27, 27, 27, 33, 33, 33, 33, 33, 33,
  33, 33, 33, 33, 33, 30, 30, 30, 30, 30, 30, 30, 33, 33, 33, 30, 30, 30, 30, 30,
  30, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 30, 30,
  30, 30, 33, 33, 33, 33, 33, 33, 30, 30, 30, 33
};

constexpr size_t twos[max_len + 1] = {
  0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2,
  2, 2, 2, 2, 3, 1, 1, 2, 2, 3, 4, 2, 2, 2, 2, 2, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10,
  10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0,
  1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0,
};

constexpr size_t code_lens[256] = {
  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
  34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
  50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
  3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11,
  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
  4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1
};

} // namespace re

class recomp_lz_helper {
  static constexpr size_t mask = 0x1fff;
  struct link {
    size_t offset;
    int next, prev;
  };
 public:
  recomp_lz_helper(std::span<const uint8_t> input)
      : input(input) {
    for (size_t i = 0; i < 256; ++i) head[i] = -1;
    for (size_t i = 0; i < 256; ++i) tail[i] = -1;
  }

  void update(std::span<const uint8_t> output, size_t i) {
    size_t j = i & mask;
    if (i > mask) {
      const uint8_t pv = output[i - mask - 1];
      int prev = link[j].prev;
      tail[pv] = prev;
      if (prev < 0) head[pv] = -1;
      else link[prev].next = -1;
    }
    const uint8_t v = output[i];
    int h = head[v];
    if (h < 0) {
      link[j] = {i, -1, -1};
      tail[v] = j;
    } else {
      link[h].prev = j;
      link[j] = {i, h, -1};
    }
    head[v] = j;
  }

  encode::lz_data find_best(
      std::span<const uint8_t> output, size_t i,
      const size_t max_ofs, const size_t max_len) const {
    const uint8_t v = input[i];
    size_t so = output.size(), si = input.size();
    encode::lz_data ret = {0, 0};
    for (int l = head[v]; l >= 0; l = link[l].next) {
      size_t o = link[l].offset, len = 1;
      if (o + max_ofs < so) break;
      assert(output[o] == input[i]);
      for (; o + len < so && i + len < si && output[o + len] == input[i + len]; ) {
        if (++len == max_len) break;
      }
      if (ret.len < len) {
        ret.ofs = o;
        ret.len = len;
        if (len == max_len) break;
      }
    }
    return ret;
  }

 private:
  std::span<const uint8_t> input;
  std::array<int, 256> head;
  std::array<int, 256> tail;
  std::array<link, mask + 1> link;
};

std::vector<uint8_t> recomp(std::span<const uint8_t> input) {
  recomp_lz_helper lz_helper(input);
  std::vector<uint8_t> ret;

  auto write8 = [&](uint8_t v) {
    ret.push_back(v);
    lz_helper.update(ret, ret.size() - 1);
  };

  auto write = [&](size_t len, size_t d) {
    if (len <= re::lzs_max_len && d <= re::lzs_max_ofs) {
      write8(0xfc | (len - 3) >> 2);
      write8((len - 3) << 6 | d);
    } else {
      write8(0xf8 | ((len - 3) >> 3));
      write8((len - 3) << 5 | d >> 8);
      write8(d);
    }
  };

  auto write8s = [&](size_t adr, size_t l) {
    for (size_t i = 0; i < l; ++i) {
      write8(input[adr + i]);
    }
  };

  size_t addr = 0;
  while (addr < input.size()) {
    uint8_t code = input[addr];
    const size_t code_len = re::code_lens[code];
    assert(code_len != 0);
    auto res_lz = lz_helper.find_best(ret, addr, re::lzl_max_ofs, re::lzl_max_len);
    size_t d = ret.size() - res_lz.ofs;
    if (res_lz.len < 3 || (res_lz.len == 3 && d > re::lzs_max_ofs)) {
      res_lz.len = 0;
    }
    if (res_lz.len == 0) {
      write8s(addr, code_len);
      addr += code_len;
    } else {
      if ((0xc0 <= code && code < 0xf0) && res_lz.len >= 3) {
        const size_t k = (d == 3) ? 3 : 0, si = input.size();
        size_t l = 3 + k;
        for (size_t i = 0; addr + i + 3 < si && input[addr + i] == input[addr + i + 3]; ++i) {
          l += 1;
          if (l == re::max_len) break;
        }
        assert(addr + l - k <= input.size());
        const size_t free_len = re::frees[l];
        assert(3 <= free_len && free_len <= l);
        if (free_len > res_lz.len) {
          write8s(addr, free_len - k);
          addr += free_len - k;
          const size_t t = re::twos[l];
          if (t > 0) {
            const size_t len2 = std::min<size_t>(free_len, 9);
            assert(t * len2 + free_len <= l);
            assert(len2 + 2 * (t - 1) <= re::lzs_max_ofs);
            for (size_t i = 0; i < t; ++i) write(len2, len2 + 2 * i);
            addr += len2 * t;
          }
          continue;
        }
      }
      write(res_lz.len, d);
      size_t l = res_lz.len;
      while (l > 0) {
        size_t clen = re::code_lens[input[addr]];
        assert(clen > 0);
        if (l >= clen) {
          l -= clen;
          addr += clen;
          continue;
        }
        write8s(addr + l, clen - l);
        addr += clen;
        break;
      }
    }
  }
  assert(addr == input.size());
  return ret;
}

} // namespace

std::vector<uint8_t> fe4_comp(std::span<const uint8_t> input) {
  enum CompType {
    uncomp, coupled,
    common_lo16, common_hi16,
    common_lo8_0, common_lo8_f, common_lo8,
    common_hi8_0, common_hi8_f, common_hi8,
    lzl, lzs,
    rle, rlel,
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());
  using cost_type = typename sssp_solver<CompType>::cost_type;

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    auto res_lz = lz_helper.find_best_closest(i, 0x3ff, 0x11);
    dp.update_lz<std::greater<cost_type>>(i, 2, 0x11, res_lz, Constant<2>(), lzs);

    auto res_lzl = lz_helper.find_best_closest(i, 0x7fff, 0x41);
    dp.update_lz<std::greater<cost_type>>(i, 0x2, 0x41, res_lzl, Constant<3>(), lzl);

    dp.update<std::greater<cost_type>>(i, 1, 0x40, Linear<1, 1>(), uncomp);

    rlen = encode::run_length(input, i, rlen);
    dp.update<std::greater<cost_type>>(i, 3, 10, rlen, Constant<2>(), rle);
    dp.update(i, 11, 0x1002, rlen, Constant<3>(), rlel);

    auto coupled8_len = encode::coupled8(input, i, 0x20);
    dp.update_k<2>(i, 2, 0x20, coupled8_len, LinearQ<1, 3, 2>(), coupled);
    auto common_lo16_len = encode::common_lo16(input, i, 0x22).len;
    dp.update_k<2>(i, 4, 0x22, common_lo16_len, LinearQ<1, 5, 2>(), common_lo16);
    auto common_hi16_len = encode::common_hi16(input, i, 0x22).len;
    dp.update_k<2>(i, 4, 0x22, common_hi16_len, LinearQ<1, 5, 2>(), common_hi16);
    auto common_hi8_res = encode::common_hi8(input, i, 0x12);
    if (((common_hi8_res.v + 0x10) & 0xf0) < 0x20) {
      dp.update(i, 3, 0x12, common_hi8_res.len,
        [](size_t i) { return (i + 4) >> 1; },
        common_hi8_res.v == 0 ? common_hi8_0 : common_hi8_f);
    } else {
      dp.update(i, 2, 0x11, common_hi8_res.len,
        [](size_t i) { return (i + 5) >> 1; }, common_hi8);
    }
    auto common_lo8_res = encode::common_lo8(input, i, 0x12);
    if (((common_lo8_res.v + 1) & 0xf) < 2) {
      dp.update(i, 3, 0x12, common_lo8_res.len,
        [](size_t i) { return (i + 4) >> 1; },
        common_lo8_res.v == 0 ? common_lo8_0 : common_lo8_f);
    } else {
      dp.update(i, 2, 0x11, common_lo8_res.len,
        [](size_t i) { return (i + 5) >> 1; }, common_lo8);
    }
    lz_helper.add_element(i);
  }
  using namespace data_type;
  writer ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x00 + cmd.len - 1, {cmd.len, &input[adr]}); break;
    case common_hi8: {
      ret.write<d8, d8, d8n_lb>(0x40 + cmd.len - 2, 0x00 | (input[adr] >> 4), {cmd.len, &input[adr]});
    } break;
    case common_hi8_0: {
      ret.write<d8, d8, d8n_lb>(0x40 + cmd.len - 3, 0x80 | (input[adr] & 0x0f), {cmd.len - 1, &input[adr + 1]});
    } break;
    case common_hi8_f: {
      ret.write<d8, d8, d8n_lb>(0x40 + cmd.len - 3, 0xc0 | (input[adr] & 0x0f), {cmd.len - 1, &input[adr + 1]});
    } break;
    case common_lo8: {
      ret.write<d8, d8, d8n_hb>(0x40 + cmd.len - 2, 0x10 | (input[adr] & 0x0f), {cmd.len, &input[adr]});
    } break;
    case common_lo8_0: {
      ret.write<d8, d8, d8n_hb>(0x40 + cmd.len - 3, 0x90 | (input[adr] >> 4), {cmd.len - 1, &input[adr + 1]});
    } break;
    case common_lo8_f: {
      ret.write<d8, d8, d8n_hb>(0x40 + cmd.len - 3, 0xd0 | (input[adr] >> 4), {cmd.len - 1, &input[adr + 1]});
    } break;
    case coupled: ret.write<d8, d8nk>(0x50 + ((cmd.len - 2) >> 1), {cmd.len, &input[adr], 2}); break;
    case common_lo16:
      ret.write<d8, d8, d8nk>(0x60 + ((cmd.len - 4) >> 1), input[adr], {cmd.len, &input[adr + 1], 2}); break;
    case common_hi16:
      ret.write<d8, d8, d8nk>(0x70 + ((cmd.len - 4) >> 1), input[adr + 1], {cmd.len, &input[adr], 2}); break;
    case lzs: ret.write<d16b>(0x8000 | ((cmd.len - 2) << 10) | (adr - cmd.lz_ofs)); break;
    case lzl: ret.write<d24b>(0xc00000 | ((cmd.len - 2) << 15) | (adr - cmd.lz_ofs)); break;
    case rlel: ret.write<d16b, d8>(0xe000 + cmd.len - 3, input[adr]); break;
    case rle: ret.write<d8, d8>(0xf0 + cmd.len - 3, input[adr]); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.total_cost() == ret.size());
  assert(adr == input.size());
  ret.write<d8>(0xff);
  ret.out = recomp(ret.out);
  return ret.out;
}

} // namespace sfc_comp
