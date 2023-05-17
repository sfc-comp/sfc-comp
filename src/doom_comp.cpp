#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> doom_comp_04(std::span<const uint8_t> in) {
  check_size(in.size(), 1, 0xffff);
  enum tag { uncomp, uncompl, lz2, lz3, lz4, lz5 };

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x104);
  auto c8 = dp.c<8>(0x108);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);

    dp.update(i, 1, 8, c8, 5, uncomp);
    dp.update(i, 9, 0x108, c8, 11, uncompl);

    dp.update(i, 2, 2, lz_helper.find(i, 0x00ff, 2), c0, 10, lz2);
    dp.update(i, 3, 3, lz_helper.find(i, 0x01ff, 3), c0, 12, lz3);
    dp.update(i, 4, 4, lz_helper.find(i, 0x03ff, 4), c0, 13, lz4);
    dp.update(i, 5, 0x104, lz_helper.find(i, 0x0fff, 5), c0, 23, lz5);

    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b16_l ret(4); ret.write<b1>(0);

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case uncomp: ret.write<bnh, bnh, b8hn>({2, 0b00}, {3, cmd.len - 1}, {cmd.len, &input[adr]}); break;
    case lz2: ret.write<bnh, bnh>({2, 0b01}, {8, d}); break;
    case lz3: ret.write<bnh, bnh>({3, 0b100}, {9, d}); break;
    case lz4: ret.write<bnh, bnh>({3, 0b101}, {10, d}); break;
    case lz5: ret.write<bnh, bnh, bnh>({3, 0b110}, {8, cmd.len - 5}, {12, d}); break;
    case uncompl: ret.write<bnh, bnh, b8hn>({3, 0b111}, {8, cmd.len - 9}, {cmd.len, &input[adr]}); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, in.size());
  write16(ret.out, 2, ret.size() - 2);
  write16(ret.out, 4, 0x8000 | (read16(ret.out, 4) >> 1));
  std::reverse(ret.out.begin() + 4, ret.out.end());
  for (size_t i = 4; i < ret.size(); i += 2) std::swap(ret.out[i], ret.out[i + 1]);
  assert(adr == in.size());
  assert(dp.optimal_cost() + 1 + 4 * 8 == ret.bit_length());
  return ret.out;
}

std::vector<uint8_t> doom_comp_08(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0xffff);
  enum tag { uncomp, rle, lzs, lz };

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x21);
  auto c1 = dp.c<1>(0x40);

  size_t rlen = 0;
  if (input.size() > 0) lz_helper.reset(input.size() - 1);
  for (size_t i = input.size(); i-- > 0; ) {
    const auto res_lz1 = lz_helper.find(i, 0x20, 1);
    if (i > 0) lz_helper.reset(i - 1);

    dp.update(i, 1, 0x40, c1, 1, uncomp);
    rlen = encode::run_length_r(input, i, rlen);
    dp.update(i, 3, 0x21, rlen, c0, 2, rle);
    dp.update(i, 1, 1, res_lz1, c0, 1, lzs);
    dp.update(i, 3, 0x12, lz_helper.find(i, 0x0800, 3), c0, 2, lz);

    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret; ret.write<d16, d8>(0, 0x01);
  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    size_t d = adr - cmd.lz_ofs();
    switch (cmd.type) {
    case rle: {
      ret.write<d8, d8>(cmd.len - 2, input[adr]);
    } break;
    case lzs: {
      ret.write<d8>(0x20 | (d - 1));
    } break;
    case uncomp: {
      ret.write<d8, d8n>(0x40 | (cmd.len - 1), {cmd.len, &in[in.size() - adr - cmd.len]});
    } break;
    case lz: {
      assert(d >= 2);
      ret.write<d16b>(0x8000 | (cmd.len - 3) << 11 | (d - 1));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, in.size());
  ret.write<d8>(0x00);
  assert(adr == in.size());
  assert(dp.optimal_cost() + 4 == ret.size());
  return ret.out;
}

namespace {

template <typename U = uint32_t>
requires std::unsigned_integral<U>
class lz_helper_doom {
public:
  using index_type = U;
  using signed_index_type = std::make_signed_t<index_type>;

public:
  lz_helper_doom(std::span<const uint8_t> in_odd, std::span<const uint8_t> in_even) {
    const size_t so = in_odd.size(), se = in_even.size();
    std::vector<int16_t> input(so + 1 + se);
    std::ranges::copy(in_odd, input.begin());
    input[so] = -1;
    std::ranges::copy(in_even, input.begin() + so + 1);

    const auto sa = suffix_array<int16_t>(input);
    const auto [lcp, rank] = sa.lcp_rank();
    this->rank = std::move(rank);
    this->lcp = decltype(this->lcp)(lcp);
    seg = decltype(seg)(rank.size());
    seg.init([&](size_t i) { return sa[i] < so ? sa[i] : seg.iden; });
  }

public:
  encode::lz_data find(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg.nodes());
  }

private:
  std::vector<index_type> rank;
  segment_tree<range_max<signed_index_type>> seg;
  segment_tree<range_min<index_type>> lcp;
};

} // namespace

std::vector<uint8_t> doom_comp_0c(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  enum tag { uncomp, rle, inc, lz, lz2 };

  std::vector<uint8_t> in_odd, in_even;
  for (size_t i = 0; i < input.size(); i += 2) in_even.push_back(input[i]);
  std::reverse(in_even.begin(), in_even.end());
  for (size_t i = 1; i < input.size(); i += 2) in_odd.push_back(input[i]);
  std::reverse(in_odd.begin(), in_odd.end());

  lz_helper helper_odd(in_odd, true), helper_even(in_even, true);
  solver<tag> dp_odd(in_odd.size()), dp_even(in_even.size());

  lz_helper_doom lz_helper_d(in_odd, in_even);

  const auto solve = [&](std::span<const uint8_t> inp,
      solver<tag>& dp, decltype(helper_odd)& lz_help, bool second) {

    auto c0 = dp.c<0>(0x22);
    auto c1 = dp.c<1>(0x40);
    size_t rlen = 0, rleni = 0;
    for (size_t i = inp.size(); i-- > 0; ) {
      lz_help.reset(i);

      dp.update(i, 1, 0x40, c1, 1, uncomp);
      rlen = encode::run_length_r(inp, i, rlen);
      dp.update(i, 3, 0x21, rlen, c0, 2, rle);
      rleni = encode::run_length_r(inp, i, rleni, 0xff);
      dp.update(i, 3, 0x22, rleni, c0, 2, inc);
      dp.update(i, 3, 0x12, lz_help.find(i, 0x800, 3), c0, 2, lz);
      if (second) {
        const auto res_lz2 = lz_helper_d.find(in_odd.size() + 1 + i, 0x0801, 3);
        dp.update(i, 3, 0x12, res_lz2, c0, 2, lz2);
      }

      c0.update(i); c1.update(i);
    }
  };
  solve(in_odd, dp_odd, helper_odd, false);
  solve(in_even, dp_even, helper_even, true);

  using namespace data_type;
  writer ret; ret.write<d16, d8>(0, 0x02);

  for (size_t k = 0; k < 2; ++k) {
    size_t adr = 0;
    const auto& dp = (k == 0) ? dp_odd : dp_even;
    const auto& in = (k == 0) ? in_odd : in_even;
    for (const auto& cmd : dp.optimal_path()) {
      switch (cmd.type) {
      case rle: {
        ret.write<d8, d8>(cmd.len - 2, in[adr]);
      } break;
      case inc: {
        ret.write<d8, d8>(0x20 | (cmd.len - 3), in[adr]);
      } break;
      case uncomp: {
        ret.write<d8>(0x40 | (cmd.len - 1));
        for (size_t i = 0; i < cmd.len; ++i) ret.write<d8>(in[adr + cmd.len - 1 - i]);
      } break;
      case lz: {
        const size_t d = adr - cmd.lz_ofs();
        ret.write<d16b>(0x8000 | (cmd.len - 3) << 11 | (d - 1));
      } break;
      case lz2: {
        const size_t d = adr + in_odd.size() - cmd.lz_ofs();
        ret.write<d16b>(0x8000 | (cmd.len - 3) << 11 | (d - 1));
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    assert(adr == in.size());
  }
  ret.write<d8>(0x00);
  write16(ret.out, 0, input.size());
  assert(3 + dp_odd.optimal_cost() + dp_even.optimal_cost() + 1 == ret.size());
  return ret.out;
}

std::vector<uint8_t> doom_comp_1(std::span<const uint8_t> input) {
  std::vector<uint8_t> ret = doom_comp_04(input);
  size_t comp_type = 1;

  if (auto res = doom_comp_08(input); res.size() < ret.size()) {
    ret = std::move(res); comp_type = 2;
  }
  if (auto res = doom_comp_0c(input); res.size() < ret.size()) {
    ret = std::move(res); comp_type = 3;
  }
  ret.resize(ret.size() + 2);
  for (size_t i = 0; i < ret.size() - 2; ++i) ret[ret.size() - 1 - i] = ret[ret.size() - 3 - i];
  write16(ret, 0, comp_type << 2);

  return ret;
}

std::vector<uint8_t> doom_comp_2(std::span<const uint8_t> input) {
  std::vector<uint8_t> ret = doom_comp_08(input);

  if (auto res = doom_comp_0c(input); res.size() < ret.size()) {
    ret = std::move(res);
  }
  return ret;
}

} // namespace sfc_comp
