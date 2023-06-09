#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

template <typename LzFunc, typename LzUpdateFunc>
requires std::convertible_to<std::invoke_result_t<LzFunc, size_t, size_t, size_t>, encode::lz_data> &&
         std::invocable<LzUpdateFunc, size_t>
std::vector<uint8_t> diet_comp_core(std::span<const uint8_t> input, const size_t header_size,
    LzFunc&& find_lz, LzUpdateFunc&& lz_helper_update) {
  enum method { uncomp, lz2, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001, 10, 0b01,       0b10},
    {0x0201, 11, 0b001,      0b100},
    {0x0401, 13, 0b00001,    0b10010},
    {0x0801, 15, 0b0000001,  0b1001010},
    {0x1001, 16, 0b00000000, 0b10010101},
  }, 0x2000);

  static constexpr auto len_tab = to_vranges({
    {0x0003,  1, 0b1},          // 1
    {0x0004,  2, 0b01},         // 01
    {0x0005,  3, 0b001},        // 001
    {0x0006,  4, 0b0001},       // 0001
    {0x0007,  6, 0b00001'0},    // 00001_
    {0x0009,  9, 0b000000'000}, // 000000___
    {0x0011, 14, 0b000001}      // 000001<1B>
  }, 0x0110);

  solver<tag> dp(input.size());
  auto c0 = dp.template c<0>(len_tab.back().max);

  if (input.size() > 0) lz_helper_update(input.size() - 1);
  for (size_t i = input.size(); i-- > 0; ) {
    dp.update(i, 1, 9, {uncomp, 0, 0});
    dp.update_matrix(i, ofs_tab, len_tab, c0, 2,
      [&](size_t oi) { return find_lz(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    if (i > 0) lz_helper_update(i - 1);
    const auto res_lz2 = find_lz(i, 0x900, 2);
    if (res_lz2.len >= 2) {
      dp.update(i, 2, 2, find_lz(i, 0x100, 2), c0, 11, {lz2, 0, 0});
      dp.update(i, 2, 2, res_lz2, c0, 14, {lz2, 1, 0});
    }
    c0.update(i);
  }

  using namespace data_type;
  writer_b16_pre_l ret(header_size);
  ret.write<none>(none());

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      ret.write<b1, d8>(true, input[adr]);
    } break;
    case lz2: {
      const size_t d = adr - cmd.lz_ofs(); assert(d >= 2);
      ret.write<bnh, d8>({2, 0}, (0x900 - d) & 0xff);
      if (oi == 0) {
        ret.write<bnh>({1, 0});
      } else {
        ret.write<bnh>({4, 0x08 | (0x900 - d) >> 8});
      }
    } break;
    case lz: {
      const auto& o = ofs_tab[oi];
      const size_t ov = (o.max - (adr - cmd.lz_ofs()));
      ret.write<bnh, d8>({2, 1}, ov & 0xff);
      ret.write<bnh>({o.bitlen - 8, masked_add(o.val, ov >> 8, o.mask)});
      const auto& l = len_tab[li];
      const size_t ld = cmd.len - l.min;
      if (l.bitlen <= 9) {
        ret.write<bnh>({l.bitlen, l.val | ld});
      } else {
        ret.write<bnh, d8>({l.bitlen - 8, l.val}, ld);
      }
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh, d8, bnh>({2, 0}, 0xff, {1, 0});
  ret.trim();

  assert(adr == input.size());
  assert(dp.optimal_cost() + 11 + header_size * 8 == ret.bit_length());

  return ret.out;
}

} // namespace

std::vector<uint8_t> diet_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  lz_helper lz_helper(input, true);
  auto ret = diet_comp_core(input, 0x11,
    [&](size_t i, size_t max_dist, size_t min_len) { return lz_helper.find(i, max_dist, min_len); },
    [&](size_t i) { lz_helper.reset(i); }
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
    [&](size_t i, size_t max_dist, size_t) { return lz_helper.find_non_overlapping(i, max_dist); },
    [&](size_t) {}
  );
  write16(ret, 0, 0x3439);
  return ret;
}

std::vector<uint8_t> pac_in_time_comp_a44a(std::span<const uint8_t> input) {
  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = to_vranges({
    {0x0002,  3, 0b100},
    {0x0003,  3, 0b101},
    {0x0004,  4, 0b00'00},
    {0x0008,  6, 0b010'000},
    {0x0010,  8, 0b0111'0000},
    {0x0020, 10, 0b01100'00000},
    {0x0040, 11, 0b01101'000000}
  }, 0x007f);
  const size_t lz_min_len = len_tab.front().min;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001 * 2 - 1,  6, 0b100'000},
    {0x0009 * 2 - 1,  8, 0b11111'000},
    {0x0011 * 2 - 1,  8, 0b1010'0000},
    {0x0021 * 2 - 1,  9, 0b1011'00000},

    {0x0041 * 2 - 1,  8, 0b00'000000},
    {0x0081 * 2 - 1, 10, 0b110'0000000},
    {0x0101 * 2 - 1, 12, 0b1110'00000000},
    {0x0201 * 2 - 1, 12, 0b011'00000000'0},
    {0x0401 * 2 - 1, 13, 0b010'00000000'00},
    {0x0801 * 2 - 1, 16, 0b11110'00000000'000}
  }, 0x0fff * 2);
  // non-decreasing parts
  static constexpr auto ofs_tab_a = std::span(ofs_tab.begin(), 4);
  static constexpr auto ofs_tab_b = std::span(ofs_tab.begin() + 4, ofs_tab.size() - 4);

  auto lz_helpers = [&] {
    lz_helper even(input); lz_helper odd(even);
    return std::to_array({std::move(even), std::move(odd)});
  }();
  std::vector<std::array<encode::lz_data, ofs_tab.size()>> lz_memo(input.size(), {{}});
  for (size_t i = 0; i < input.size(); ++i) {
    const auto f = [&](std::span<const vrange> o_tab, size_t beg) {
      if (size_t j = i + o_tab.front().min; j < input.size()) {
        encode::lz::find_all(j, o_tab, lz_min_len, std::span(lz_memo[j].data() + beg, o_tab.size()),
          [&](size_t max_ofs) { return lz_helpers[j & 1].find(j, max_ofs, lz_min_len); }
        );
      }
    };
    lz_helpers[i & 1].add_element(i);
    f(ofs_tab_a, 0);
    f(ofs_tab_b, ofs_tab_a.size());
  }

  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(len_tab.back().max);
  for (size_t i = input.size(); i-- > 0; ) {
    dp.update(i, 1, 10, {uncomp, 0, 0});
    const auto f = [&](std::span<const vrange> o_tab, size_t beg) {
      dp.update_matrix(i, o_tab, len_tab, c0, 0,
        [&](size_t oi) { return lz_memo[i][oi + beg]; },
        [&](size_t oi, size_t li) -> tag { return {lz, oi + beg, li}; }
      );
    };
    f(ofs_tab_a, 0);
    f(ofs_tab_b, ofs_tab_a.size());
    c0.update(i);
  }

  using namespace data_type;
  writer_b16_pre_l ret(2); ret.write<none>(none());

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      ret.write<bnh, d8>({2, 0b11}, input[adr]);
    } break;

    case lz: {
      const size_t d = adr - cmd.lz_ofs(); assert(d % 2 == 0);
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[oi];
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
  assert(dp.optimal_cost() + 2 * 8 +
         len_tab.front().bitlen + ofs_tab.back().bitlen == ret.bit_length());
  write16(ret.out, 0, 0xa44a);
  return ret.out;
}

std::vector<uint8_t> pac_in_time_comp_a55a(std::span<const uint8_t> input) {
  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = to_vranges({
    {0x0002,  1, 0b0},
    {0x0003,  4, 0b1111},
    {0x0004,  5, 0b100'00},
    {0x0008,  6, 0b101'000},
    {0x0010,  7, 0b110'0000},
    {0x0020, 10, 0b11100'00000},
    {0x0040, 11, 0b11101'000000}
  }, 0x007f);

  static constexpr auto ofs_tab = to_vranges({
    {0x0002,  8, 0b0000000'0},
    {0x0081, 10, 0b0000000'1'00},
    {0x0101, 11, 0b0000000'1'01'0},
    {0x0201, 13, 0b0000000'1'100'00},
    {0x0401, 14, 0b0000000'1'101'000},
    {0x0801, 15, 0b0000000'1'110'0000},
    {0x1001, 16, 0b0000000'1'111'00000},
  }, 0x2000);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.c<0>(len_tab.back().max);

  if (input.size() > 0) lz_helper.reset(input.size() - 1);
  for (size_t i = input.size(); i-- > 0; ) {
    if (i > 0) lz_helper.reset(i - 1);
    dp.update(i, 1, 9, {uncomp, 0, 0});
    dp.update_matrix(i, ofs_tab, len_tab, c0, 1,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i);
  }

  using namespace data_type;
  writer_b16_pre_l ret(2); ret.write<none>(none());

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      ret.write<b1, d8>(true, input[adr]);
    } break;

    case lz: {
      const size_t d = adr - cmd.lz_ofs(); assert(d >= 2);
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[oi];
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
  assert(dp.optimal_cost() + 9 + 2 * 8 == ret.bit_length());

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
