#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

template <typename LzFunc, typename LzUpdateFunc>
requires std::convertible_to<std::invoke_result_t<LzFunc, size_t, size_t>, encode::lz_data> &&
         std::invocable<LzUpdateFunc, size_t>
std::vector<uint8_t> diet_comp_core(std::span<const uint8_t> input, const size_t header_size,
    LzFunc&& find_lz, LzUpdateFunc&& lz_helper_update) {
  enum method { uncomp, lz2, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = std::to_array<vrange>({
    vrange(0x0001, 0x0200, 10, 0b01),       // _1
    vrange(0x0201, 0x0400, 11, 0b001),      // _01
    vrange(0x0401, 0x0800, 13, 0b00001),    // _00_1
    vrange(0x0801, 0x1000, 15, 0b0000001),  // _00_0_1
    vrange(0x1001, 0x2000, 16, 0b00000000), // _00_0_0_
  });

  static constexpr auto len_tab = std::to_array<vrange>({
    vrange(0x0003, 0x0003,  1, 0b1),          // 1
    vrange(0x0004, 0x0004,  2, 0b01),         // 01
    vrange(0x0005, 0x0005,  3, 0b001),        // 001
    vrange(0x0006, 0x0006,  4, 0b0001),       // 0001
    vrange(0x0007, 0x0008,  6, 0b00001'0),    // 00001_
    vrange(0x0009, 0x0010,  9, 0b000000'000), // 000000___
    vrange(0x0011, 0x0110, 14, 0b000001)      // 000001<1B>
  });

  sssp_solver<tag> dp(input.size());

  encode::lz_data res_lz2 = {}, res_lz2s = {};
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    if (res_lz2.len >= 2) {
      if (res_lz2s.len >= 2) {
        dp.update_lz(i, 2, 2, res_lz2s, Constant<11>(), {lz2, 0, 0});
      } else {
        dp.update_lz(i, 2, 2, res_lz2, Constant<14>(), {lz2, 1, 0});
      }
    }
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return find_lz(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      2
    );
    // dist should be >= 2.
    if (i + 1 < input.size()) {
      res_lz2 = find_lz(i + 1, 0x900);
      if (res_lz2.len >= 2) {
        res_lz2s = find_lz(i + 1, 0x100);
      }
    }
    lz_helper_update(i);
  }

  using namespace data_type;
  writer_b16_pre_l ret(header_size);
  ret.write<none>(none());

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, d8>(true, input[adr]);
    } break;
    case lz2: {
      const size_t d = adr - cmd.lz_ofs; assert(d >= 2);
      ret.write<bnh, d8>({2, 0}, (0x900 - d) & 0xff);
      if (cmd.type.oi == 0) {
        ret.write<bnh>({1, 0});
      } else {
        ret.write<bnh>({4, 0x08 | (0x900 - d) >> 8});
      }
    } break;
    case lz: {
      const size_t d = adr - cmd.lz_ofs;
      const size_t oi = cmd.type.oi;
      const auto& o = ofs_tab[oi];

      ret.write<bnh, d8>({2, 1}, -d & 0xff);
      size_t v = (o.max - d) >> 8;
      if (oi == 0) {
        v <<= 1;
      } else if (oi == 1) {
        v <<= 2;
      } else if (oi == 2) {
        v <<= 1;
        v += (v & 0x04) * 3;
      } else if (oi == 3) {
        v <<= 1;
        v += v & ~3;
        v += (v & 0x10) * 3;
      } else {
        v += v & ~1;
        v += v & ~7;
        v += (v & 0x20) * 3;
      }
      ret.write<bnh>({o.bitlen - 8, v | o.val});

      const auto& l = len_tab[cmd.type.li];
      const size_t len_v = cmd.len - l.min;
      if (l.bitlen < 10) {
        ret.write<bnh>({l.bitlen, len_v | l.val});
      } else {
        ret.write<bnh, d8>({l.bitlen - 8, l.val}, len_v);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh, d8, bnh>({2, 0}, 0xff, {1, 0});
  ret.trim();

  assert(adr == input.size());
  assert(dp.total_cost() + 11 + header_size * 8 == ret.bit_length());

  return ret.out;
}

} // namespace

std::vector<uint8_t> diet_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  lz_helper lz_helper(input);
  auto ret = diet_comp_core(input, 0x11,
    [&](size_t i, size_t max_dist) { return lz_helper.find_best(i, max_dist); },
    [&](size_t i) { lz_helper.add_element(i); }
  );
  static constexpr auto header = std::to_array<uint8_t>({
    0xb4, 0x4c, 0xcd, 0x21, 0x9d, 0x89, 0x64, 0x6c, 0x7a
  });
  std::ranges::copy(header, ret.begin());
  ret[9] = ((ret.size() - 0x11) >> 16) & 0x0f;
  write16(ret, 10, ret.size() - 0x11);
  write16(ret, 12, utility::crc16(ret, 0x11, ret.size() - 0x11));
  ret[14] = ((input.size() >> 16) & 0x3f) << 2;
  write16(ret, 15, input.size());
  return ret;
}

namespace {

std::vector<uint8_t> pac_in_time_comp_3439(std::span<const uint8_t> input) {
  non_overlapping_lz_helper lz_helper(input);
  auto ret = diet_comp_core(input, 2,
    [&](size_t i, size_t max_dist) { return lz_helper.find_non_overlapping(i, max_dist); },
    [&](size_t) {}
  );
  write16(ret, 0, 0x3439);
  return ret;
}

std::vector<uint8_t> pac_in_time_comp_a44a(std::span<const uint8_t> input) {
  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = std::to_array<vrange>({
    vrange(0x0002, 0x0002,  3, 0b100),
    vrange(0x0003, 0x0003,  3, 0b101),
    vrange(0x0004, 0x0007,  4, 0b00'00),
    vrange(0x0008, 0x000f,  6, 0b010'000),
    vrange(0x0010, 0x001f,  8, 0b0111'0000),
    vrange(0x0020, 0x003f, 10, 0b01100'00000),
    vrange(0x0040, 0x007f, 11, 0b01101'000000)
  });
  const size_t lz_min_len = len_tab.front().min;

  static constexpr auto ofs_tab = std::to_array<vrange>({
    vrange(0x0001 * 2 - 1, 0x0008 * 2,  6, 0b100'000),
    vrange(0x0009 * 2 - 1, 0x0010 * 2,  8, 0b11111'000),
    vrange(0x0011 * 2 - 1, 0x0020 * 2,  8, 0b1010'0000),
    vrange(0x0021 * 2 - 1, 0x0040 * 2,  9, 0b1011'00000),

    vrange(0x0041 * 2 - 1, 0x0080 * 2,  8, 0b00'000000),
    vrange(0x0081 * 2 - 1, 0x0100 * 2, 10, 0b110'0000000),
    vrange(0x0101 * 2 - 1, 0x0200 * 2, 12, 0b1110'00000000),
    vrange(0x0201 * 2 - 1, 0x0400 * 2, 12, 0b011'00000000'0),
    vrange(0x0401 * 2 - 1, 0x0800 * 2, 13, 0b010'00000000'00),
    vrange(0x0801 * 2 - 1, 0x0fff * 2, 16, 0b11110'00000000'000)
  });
  // non-decreasing parts
  static constexpr auto ofs_tab_a = std::span(ofs_tab.begin(), 4);
  static constexpr auto ofs_tab_b = std::span(ofs_tab.begin() + 4, ofs_tab.size() - 4);

  auto lz_helpers = [&] {
    lz_helper even(input); lz_helper odd(even);
    return std::to_array({std::move(even), std::move(odd)});
  }();
  sssp_solver<tag> dp(input.size());

  const size_t memo_size = std::bit_ceil(ofs_tab_b.front().min);
  std::array<std::array<encode::lz_data, ofs_tab_b.size()>, memo_size> lz_memo = {};

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<10>(), {uncomp, 0, 0});
    dp.update_lz_matrix(i, ofs_tab_a, len_tab,
      [&](size_t oi) { return lz_helpers[i & 1].find_best(i, ofs_tab_a[oi].max); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      0
    );
    dp.update_lz_matrix(i, ofs_tab_b, len_tab,
      [&](size_t oi) { return lz_memo[i % memo_size][oi]; },
      [&](size_t oi, size_t li) -> tag { return {lz, oi + ofs_tab_a.size(), li}; },
      0
    );
    lz_helpers[i & 1].add_element(i);
    if (size_t j = i + ofs_tab_b.front().min; j < input.size()) {
      lz::find_all(j, ofs_tab_b, lz_min_len, lz_memo[j % memo_size],
        [&](size_t max_ofs) { return lz_helpers[j & 1].find_best(j, max_ofs); }
      );
    }
  }

  using namespace data_type;
  writer_b16_pre_l ret(2); ret.write<none>(none());

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<bnh, d8>({2, 0b11}, input[adr]);
    } break;

    case lz: {
      const size_t d = adr - cmd.lz_ofs; assert(d % 2 == 0);
      const auto& l = len_tab[cmd.type.li];
      const auto& o = ofs_tab[cmd.type.oi];
      assert(l.min <= cmd.len && cmd.len <= l.max);
      assert(o.min <= d && d <= o.max);
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      size_t v = (std::bit_ceil(o.max) - d) / 2;
      if ((o.max - o.min) / 2 + 1 >= 0x100) {
        size_t bits = ilog2(d / 2 - 1); assert(bits >= 8);
        size_t prefix_bits = o.bitlen - bits;
        ret.write<bnh>({prefix_bits, o.val >> bits});
        ret.write<d8>(v >> (bits - 8));
        ret.write<bnh>({bits - 8, v & low_bits_mask(bits - 8)});
      } else {
        ret.write<bnh>({o.bitlen, o.val + v});
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({len_tab.front().bitlen, len_tab.front().val});
  ret.write<bnh>({5, ofs_tab.back().val >> 11});
  ret.write<d8>(0);
  ret.write<bnh>({3, 0});
  ret.trim();

  assert(adr == input.size());
  assert(dp.total_cost() + 2 * 8 +
         len_tab.front().bitlen + ofs_tab.back().bitlen == ret.bit_length());
  write16(ret.out, 0, 0xa44a);
  return ret.out;
}

std::vector<uint8_t> pac_in_time_comp_a55a(std::span<const uint8_t> input) {
  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = std::to_array<vrange>({
    vrange(0x0002, 0x0002,  1, 0b0),
    vrange(0x0003, 0x0003,  4, 0b1111),
    vrange(0x0004, 0x0007,  5, 0b100'00),
    vrange(0x0008, 0x000f,  6, 0b101'000),
    vrange(0x0010, 0x001f,  7, 0b110'0000),
    vrange(0x0020, 0x003f, 10, 0b11100'00000),
    vrange(0x0040, 0x007f, 11, 0b11101'000000)
  });

  static constexpr auto ofs_tab = std::to_array<vrange>({
    vrange(0x0002, 0x0080,  8, 0b0000000'0),
    vrange(0x0081, 0x0100, 10, 0b0000000'1'00),
    vrange(0x0101, 0x0200, 11, 0b0000000'1'01'0),
    vrange(0x0201, 0x0400, 13, 0b0000000'1'100'00),
    vrange(0x0401, 0x0800, 14, 0b0000000'1'101'000),
    vrange(0x0801, 0x1000, 15, 0b0000000'1'110'0000),
    vrange(0x1001, 0x2000, 16, 0b0000000'1'111'00000),
  });

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      1
    );
    if (i > 0) lz_helper.add_element(i - 1);
  }

  using namespace data_type;
  writer_b16_pre_l ret(2); ret.write<none>(none());

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, d8>(true, input[adr]);
    } break;

    case lz: {
      const size_t d = adr - cmd.lz_ofs; assert(d >= 2);
      const auto& l = len_tab[cmd.type.li];
      const auto& o = ofs_tab[cmd.type.oi];
      ret.write<b1, bnh>(false, {l.bitlen, l.val + (cmd.len - l.min)});
      size_t v = o.max - d;
      if (d <= ofs_tab[0].max) {
        ret.write<d8>(v << 1);
      } else {
        size_t bits = ilog2(d - 1); assert(bits >= 7);
        ret.write<d8>((v >> (bits - 7)) << 1 | 1);
        ret.write<bnh>({o.bitlen - 8, o.val + (v & low_bits_mask(bits - 7))});
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b1, d8>(false, 0xfe);
  ret.trim();
  assert(adr == input.size());
  assert(dp.total_cost() + 9 + 2 * 8 == ret.bit_length());

  write16(ret.out, 0, 0xa55a);
  return ret.out;
}

} // namespace

std::vector<uint8_t> pac_in_time_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);
  std::vector<uint8_t> best;
  for (const auto& comp : {pac_in_time_comp_a55a, pac_in_time_comp_a44a, pac_in_time_comp_3439}) {
    if (auto res = comp(input); best.empty() || res.size() < best.size()) {
      best = std::move(res);
    }
  }
  return best;
}

} // namespace sfc_comp
