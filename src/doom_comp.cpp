#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> doom_comp_04(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);
  enum CompType {
    uncomp, uncompl,
    lz2, lz3, lz4, lz5
  };

  std::vector<uint8_t> input2(input.rbegin(), input.rend());

  lz_helper lz_helper(input2);
  sssp_solver<CompType> dp(input2.size());

  for (size_t i = 0; i < input2.size(); ++i) {
    dp.update(i, 1, 8, Linear<8, 5>(), uncomp);
    dp.update(i, 9, 0x108, Linear<8, 11>(), uncompl);

    auto res_lz2 = lz_helper.find_best(i, 0x00ff);
    dp.update_lz(i, 2, 2, res_lz2, Constant<10>(), lz2);
    auto res_lz3 = lz_helper.find_best(i, 0x01ff);
    dp.update_lz(i, 3, 3, res_lz3, Constant<12>(), lz3);
    auto res_lz4 = lz_helper.find_best(i, 0x03ff);
    dp.update_lz(i, 4, 4, res_lz4, Constant<13>(), lz4);
    auto res_lz5 = lz_helper.find_best(i, 0x0fff);
    dp.update_lz(i, 5, 0x104, res_lz5, Constant<23>(), lz5);

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16 ret; ret.write<d16, d16>(0, 0);

  ret.write<b1l>(0);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: {
      ret.write<b8ln_h>({5, cmd.len - 1});
      for (size_t i = 0; i < cmd.len; ++i) ret.write<b8ln_h>({8, input2[adr + i]});
    } break;
    case lz2: {
      ret.write<b8ln_h, b8ln_h>({2, 1}, {8, d});
    } break;
    case lz3: {
      ret.write<b8ln_h, b8ln_h>({3, 4}, {9, d});
    } break;
    case lz4: {
      ret.write<b8ln_h, b8ln_h>({3, 5}, {10, d});
    } break;
    case lz5: {
      ret.write<b8ln_h, b8ln_h, b8ln_h>({3, 6}, {8, cmd.len - 5}, {12, d});
    } break;
    case uncompl: {
      ret.write<b8ln_h>({11, 0x700 | (cmd.len - 9)});
      for (size_t i = 0; i < cmd.len; ++i) ret.write<b8ln_h>({8, input2[adr + i]});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  write16(ret.out, 2, ret.size() - 2);
  write16(ret.out, 4, 0x8000 | (read16(ret.out, 4) >> 1));
  std::reverse(ret.out.begin() + 4, ret.out.end());
  for (size_t i = 4; i < ret.size(); i += 2) std::swap(ret.out[i], ret.out[i + 1]);
  assert(dp.total_cost() + 1 + 4 * 8 == ret.bit_length());
  assert(adr == input.size());
  return ret.out;
}

std::vector<uint8_t> doom_comp_08(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  enum CompType {
    uncomp,
    rle, lzs, lz
  };

  std::vector<uint8_t> input2(input.rbegin(), input.rend());

  lz_helper lz_helper(input2);
  sssp_solver<CompType> dp(input2.size());

  size_t rlen = 0;
  encode::lz_data curr_lz = {0, 0};
  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 0x40, Linear<1, 1>(), uncomp);
    rlen = encode::run_length(input2, i, rlen);
    dp.update(i, 3, 0x21, rlen, Constant<2>(), rle);

    auto res_lzs = lz_helper.find_best(i, 0x20);
    dp.update_lz(i, 1, 1, res_lzs, Constant<1>(), lzs);
    dp.update_lz(i, 3, 0x12, curr_lz, Constant<2>(), lz);
    if (i + 1 < input.size()) {
      curr_lz = lz_helper.find_best(i + 1, 0x0800);
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer ret; ret.write<d16, d8>(0, 0x01);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case rle: {
      ret.write<d8, d8>(cmd.len - 2, input2[adr]);
    } break;
    case lzs: {
      ret.write<d8>(0x20 | (d - 1));
    } break;
    case uncomp: {
      ret.write<d8, d8n>(0x40 | (cmd.len - 1), {cmd.len, &input[input.size() - adr - cmd.len]});
    } break;
    case lz: {
      assert(d >= 2);
      ret.write<d16b>(0x8000 | (cmd.len - 3) << 11 | (d - 1));
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  write16(ret.out, 0, input.size());
  ret.write<d8>(0x00);
  assert(dp.total_cost() + 4 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

namespace {

class lz_helper_doom {
public:
  lz_helper_doom(std::span<const uint8_t> in_odd, std::span<const uint8_t> in_even) : n(in_odd.size()) {
    const size_t so = in_odd.size(), se = in_even.size();
    std::vector<int16_t> input(so + 1 + se);
    for (size_t i = 0; i < so; ++i) input[i] = in_odd[i];
    input[so] = -1;
    for (size_t i = 0; i < se; ++i) input[i + so + 1] = in_even[i];

    const auto plcp = suffix_array<int16_t>(input).compute_lcp_and_rank(input);
    seg_lcp = segment_tree<range_min<int>>(plcp.first);
    rank = std::move(plcp.second);
    seg = segment_tree<range_max<int>>(input.size());

    for (size_t i = 0; i < so; ++i) seg.update(rank[i], i);
  }

public:
  encode::lz_data find_best(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank[pos], block_size, seg_lcp, seg);
  }

private:
  const size_t n;
  std::vector<int> rank;
  segment_tree<range_max<int>> seg;
  segment_tree<range_min<int>> seg_lcp;
};

} // namespace

std::vector<uint8_t> doom_comp_0c(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  enum CompType {
    uncomp,
    rle, inc, lz, lz2
  };

  std::vector<uint8_t> in_odd, in_even;
  for (size_t i = 0; i < input.size(); i += 2) in_even.push_back(input[i]);
  for (size_t i = 1; i < input.size(); i += 2) in_odd.push_back(input[i]);
  std::reverse(in_even.begin(), in_even.end());
  std::reverse(in_odd.begin(), in_odd.end());

  lz_helper helper_odd(in_odd), helper_even(in_even);
  sssp_solver<CompType> dp_odd(in_odd.size()), dp_even(in_even.size());

  lz_helper_doom lz_helper_d(in_odd, in_even);

  auto solve = [&](std::span<const uint8_t> inp, sssp_solver<CompType>& dp, lz_helper& lz_help, bool is_second) {
    size_t rlen = 0, rleni = 0;
    for (size_t i = 0; i < inp.size(); ++i) {
      dp.update(i, 1, 0x40, Linear<1, 1>(), uncomp);

      rlen = encode::run_length(inp, i, rlen);
      dp.update(i, 3, 0x21, rlen, Constant<2>(), rle);

      rleni = encode::run_length(inp, i, rleni, 0xff);
      dp.update(i, 3, 0x22, rleni, Constant<2>(), inc);

      auto res_lz = lz_help.find_best(i, 0x800);
      dp.update_lz(i, 3, 0x12, res_lz, Constant<2>(), lz);

      if (is_second) {
        auto res_lz2 = lz_helper_d.find_best(in_odd.size() + 1 + i, 0x0801);
        dp.update_lz(i, 3, 0x12, res_lz2, Constant<2>(), lz2);
      }

      lz_help.add_element(i);
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
    for (const auto cmd : dp.commands()) {
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
        size_t d = adr - cmd.lz_ofs;
        ret.write<d16b>(0x8000 | (cmd.len - 3) << 11 | (d - 1));
      } break;
      case lz2: {
        size_t d = adr + in_odd.size() - cmd.lz_ofs;
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
  assert(3 + dp_odd.total_cost() + dp_even.total_cost() + 1 == ret.size());
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
