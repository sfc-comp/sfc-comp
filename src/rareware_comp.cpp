#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

struct pre_table {
  pre_table(std::span<const uint8_t> input) {
    // heuristic
    std::vector<size_t> b8_counts(256);
    for (size_t i = 0; i < input.size(); ++i) {
      b8_counts[input[i]] += 1;
    }
    b8 = utility::k_most<uint8_t, 2>(b8_counts);
    std::fill_n(b8_counts.data(), 256, 0);
    for (size_t i = 0; i < input.size(); ) {
      size_t rlen = encode::run_length(input, i, 0);
      if (rlen >= 3) b8_counts[input[i]] += 1;
      i += rlen;
    }
    rle_b8 = utility::k_most<uint8_t, 2>(b8_counts);
  }
  std::array<uint8_t, 2> b8 = {};
  std::array<uint8_t, 2> rle_b8 = {};
};

} // namespace

std::vector<uint8_t> rareware_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum CompType {
    uncomp, uncomp1, uncomp2, rle,
    pre8_1, pre8_2, pre_rle8_1, pre_rle8_2,
    lzvs, lzs, lzm, lzl,
    pre16_1, pre16s,
    prev8, prev16
  };

  const size_t num_candidates[] = {
    1024, 512, 256, 128, 96, 64, 48, 32, 25, 20, 17};
  const size_t phase_total = sizeof(num_candidates) / sizeof(*num_candidates);

  std::vector<encode::lz_data> lzl_memo(input.size());
  std::vector<encode::lz_data> lzm_memo(input.size());
  std::vector<std::array<encode::lz_data, 0x10>> lz_memo(input.size());
  {
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      lzl_memo[i] = lz_helper.find_best(i, 0xffff);
      const size_t j = i + 0x0102;
      if (j < input.size()) lzm_memo[j] = lz_helper.find_best(j, 0x1102);
      for (size_t l = 3; l <= 0x12; ++l) {
        const size_t j = i + l - 1;
        if (j >= input.size()) break;
        lz_memo[j][l - 3] = lz_helper.find_best(j, l + 0xff);
      }
      lz_helper.add_element(i);
    }
  }

  auto candidate = utility::k_most_freq_u16(input, num_candidates[0]);
  std::vector<int64_t> pre16(0x10000, -1);
  pre_table pre(input);

  for (size_t phase = 0; phase < phase_total; ++phase) {
    for (size_t i = 0; i < candidate.size(); ++i) pre16[candidate[i]] = i;

    sssp_solver<CompType> dp(input.size());

    size_t rlen = 0;
    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, Constant<3>(), uncomp1);
      dp.update(i, 2, 2, Constant<5>(), uncomp2);
      dp.update(i, 3, 0x0f, Linear<2, 2>(), uncomp);

      rlen = encode::run_length(input, i, rlen);
      if (input[i] == pre.rle_b8[0]) {
        dp.update(i, 3, 0x12, rlen, Constant<2>(), pre_rle8_1);
      } else if (input[i] == pre.rle_b8[1]) {
        dp.update(i, 3, 0x12, rlen, Constant<2>(), pre_rle8_2);
      } else {
        dp.update(i, 3, 0x12, rlen, Constant<4>(), rle);
      }

      if (input[i] == pre.b8[0]) {
        dp.update(i, 1, 1, Constant<1>(), pre8_1);
      } else if (input[i] == pre.b8[1]) {
        dp.update(i, 1, 1, Constant<1>(), pre8_2);
      }

      if (i >= 1 && input[i] == input[i - 1]) {
        dp.update(i, 1, 1, Constant<1>(), prev8);
      }
      if (i + 1 < input.size()) {
        const uint16_t v16 = read16(input, i);
        if (i >= 2 && read16(input, i - 2) == v16) {
          dp.update(i, 2, 2, Constant<1>(), prev16);
        }
        for (size_t ofs = 2; ofs <= std::min<size_t>(i, 0x11); ++ofs) {
          if (read16(input, i - ofs) == v16) {
            dp.update_lz(i, 2, 2, encode::lz_data(i - ofs, 2), Constant<2>(), lzvs);
            break;
          }
        }
        const int16_t ind = pre16[v16];
        if (ind == 0) dp.update(i, 2, 2, Constant<1>(), pre16_1);
        else if (ind > 0) dp.update_lz(i, 2, 2, encode::lz_data(ind, 2), Constant<2>(), pre16s);
      }
      for (size_t l = 3; l <= 0x12; ++l) {
        dp.update_lz(i, l, l, lz_memo[i][l - 3], Constant<4>(), lzs);
      }
      dp.update_lz(i, 3, 0x12, lzm_memo[i], Constant<5>(), lzm);
      dp.update_lz(i, 3, 0x12, lzl_memo[i], Constant<6>(), lzl);
    }

    if (phase + 1 < phase_total) {
      for (size_t i = 0; i < candidate.size(); ++i) pre16[candidate[i]] = 0;
      for (const auto cmd : dp.commands()) {
        if (cmd.type == pre16_1) pre16[candidate[0]] += 1;
        else if (cmd.type == pre16s) pre16[candidate[cmd.lz_ofs]] += 2;
      }
      const size_t next_k = num_candidates[phase + 1];
      std::partial_sort(
        candidate.begin(), candidate.begin() + next_k, candidate.end(),
        [&](const uint16_t a, const uint16_t b) { return pre16[a] > pre16[b]; });
      for (size_t i = next_k; i < candidate.size(); ++i) pre16[candidate[i]] = -1;
      candidate.resize(next_k);
    } else {
      using namespace data_type;
      writer_b4 ret;
      for (size_t i = 0; i < 0x27; ++i) ret.write<d8>(0);
      size_t adr = 0;
      for (const auto cmd : dp.commands()) {
        size_t d = adr - cmd.lz_ofs;
        switch (cmd.type) {
        case uncomp: {
          ret.write<d8, d8n>(cmd.len, {cmd.len, &input[adr]});
        } break;
        case uncomp1: {
          ret.write<b4h, d8>(1, input[adr]);
        } break;
        case uncomp2: {
          ret.write<b4h, d8, d8>(2, input[adr], input[adr + 1]);
        } break;
        case rle: {
          ret.write<b4h, b4h, d8>(3, cmd.len - 3, input[adr]);
        } break;
        case pre_rle8_1: {
          ret.write<b4h, b4h>(4, cmd.len - 3);
        } break;
        case pre_rle8_2: {
          ret.write<b4h, b4h>(5, cmd.len - 3);
        } break;
        case pre16_1: {
          ret.write<b4h>(6);
        } break;
        case pre8_1: {
          ret.write<b4h>(7);
        } break;
        case pre8_2: {
          ret.write<b4h>(8);
        } break;
        case lzvs: {
          assert(d >= 2);
          ret.write<b4h, b4h>(9, d - 2);
        } break;
        case lzs: {
          assert(d >= cmd.len);
          assert(d - cmd.len < 0x100);
          ret.write<b4h, b4h, d8>(10, cmd.len - 3, d - cmd.len);
        } break;
        case lzm: {
          assert(d >= 0x103 && d < 0x1103);
          ret.write<b4h, b4h, b4h, d8>(11, cmd.len - 3, (d - 0x103) >> 8, d - 0x103);
        } break;
        case lzl: {
          ret.write<b4h, b4h, d16b>(12, cmd.len - 3, d);
        } break;
        case prev8: {
          ret.write<b4h>(13);
        } break;
        case prev16: {
          ret.write<b4h>(14);
        } break;
        case pre16s: {
          ret.write<b4h, b4h>(15, cmd.lz_ofs - 1);
        } break;
        }
        adr += cmd.len;
      }
      assert(0x27 + (dp.total_cost() + 1) / 2 == ret.size());
      assert(adr == input.size());
      ret.write<d8>(0x00);
      for (size_t i = 0; i < 2; ++i) ret[i + 1] = pre.rle_b8[i];
      for (size_t i = 0; i < 2; ++i) ret[i + 3] = pre.b8[i];
      for (size_t i = 0; i < 17; ++i) write16(ret.out, 2 * i + 5, candidate[i]);
      return ret.out;
    }
  }

  throw std::logic_error("phase_total == 0");
}

} // namespace sfc_comp
