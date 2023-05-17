#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> chrono_trigger_comp_fast_core(std::span<const uint8_t> input, const size_t max_len_bits) {
  check_size(input.size(), 0, 0x10000);

  enum tag { uncomp, lz };

  std::vector<uint8_t> best;

  for (size_t len_bits = 4; len_bits <= max_len_bits; ++len_bits) {
    const size_t lz_min_len = 3;
    const size_t lz_max_len = ((1 << len_bits) - 1) + lz_min_len;
    const size_t lz_max_ofs = (0x10000 >> len_bits) - 1;

    lz_helper lz_helper(input, true);
    solver<tag> dp(input.size());
    auto c0 = dp.c<0>(lz_max_len);

    for (size_t i = input.size(); i-- > 0; ) {
      lz_helper.reset(i);
      dp.update(i, 1, 9, uncomp);
      dp.update(i, lz_min_len, lz_max_len,
                lz_helper.find(i, lz_max_ofs, lz_min_len), c0, 17, lz);
      c0.update(i);
    }

    using namespace data_type;
    writer_b8_l ret(2);
    size_t adr = 0;
    for (const auto& cmd : dp.optimal_path()) {
      switch (cmd.type) {
      case uncomp: ret.write<b1, d8>(false, input[adr]); break;
      case lz: ret.write<b1, d16>(true, (adr - cmd.lz_ofs()) | (cmd.len - lz_min_len) << (16 - len_bits)); break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    assert(adr == input.size());
    assert(dp.optimal_cost() + 2 * 8 == ret.bit_length());

    const size_t method_bit = (len_bits == 4) ? 0x00 : 0x40;
    if (ret.bit == 0) {
      write16(ret.out, 0, ret.size() - 2);
    } else {
      const size_t bits_pos = ret.bits_pos;
      write16(ret.out, 0, bits_pos - 2);

      const size_t len = (8 - ret.bit) + std::popcount(ret.out[bits_pos]) + 1;
      assert(ret.size() == bits_pos + len);
      ret.write<d8, d8, d8>(0, 0, 0);

      for (size_t i = 0; i < len; ++i) ret[ret.size() - 1 - i] = ret[ret.size() - 4 - i];
      ret[bits_pos] = (8 - ret.bit) | method_bit;
      write16(ret.out, bits_pos + 1, ret.size());
      ret[bits_pos + 3] |= low_bits_mask(ret.bit) << (8 - ret.bit); // Avoids 0x00. (cf. $C3:07B7, $C3:0879, etc. in Chrono Trigger)
    }
    ret.write<d8>(method_bit);
    if (best.empty() || ret.size() < best.size()) best = std::move(ret.out);
  }
  return best;
}

std::vector<uint8_t> chrono_trigger_comp_core(
    std::span<const uint8_t> input, const size_t max_len_bits, const size_t max_bits) {
  check_size(input.size(), 0, 0x10000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  std::vector<uint8_t> best;

  for (size_t len_bits = 4; len_bits <= max_len_bits; ++len_bits) {
    const size_t lz_min_len = 3;
    const size_t lz_max_len = ((1 << len_bits) - 1) + lz_min_len;
    const size_t lz_max_ofs = (0x10000 >> len_bits) - 1;

    lz_helper lz_helper(input, true);
    std::array<solver<tag>, 8> dp;
    for (size_t b = 0; b < 8; ++b) {
      dp[b] = solver<tag>(input.size(), (b == 0) ? input.size() : -1);
    }
    auto c0s = create_array<decltype(dp[0].c<0>(0)), 8>([&](size_t b) {
      return dp[b].c<0>(lz_max_len);
    });
    auto c1 = dp[0].c<1>(lz_max_len + max_bits);

    for (size_t i = input.size(); i-- > 0; ) {
      lz_helper.reset(i);
      const auto res_lz = lz_helper.find(i, lz_max_ofs, lz_min_len);
      const size_t lz_len = std::min(res_lz.len, lz_max_len);
      const auto update = [&](size_t b, size_t to, size_t c) {
        dp[b].update_c(i, 1, c0s[to][i + 1] + c + 1, {uncomp, to, 0});
        dp[b].update(i, lz_min_len, lz_max_len, res_lz, c0s[to], c + 2, {lz, to, 0});
      };
      for (size_t b = 0; b < 8; ++b) {
        const size_t c = (b == 0) ? 1 : 0;
        update(b, (b - 1) & 7, c);
        if (b != 1) update(b, 0, c + 3);
        if (res_lz.len >= lz_min_len) {
          const size_t remain = max_bits - 1 - ((8 - b) & 7);
          const auto e = c1.find(i, lz_len + 1, lz_len + remain);
          if (e.len == c1.nlen) continue;
          const size_t u = e.len - lz_len;
          dp[b].update_c(i, e.len, (3 + c) + 2 + (e.cost - lz_len), {lz, 0, u}, res_lz.ofs);
        }
      }
      for (size_t b = 0; b < 8; ++b) c0s[b].update(i);
      c1.update(i);
    }

    const size_t method_bit = (len_bits == 4) ? 0x00 : 0x40;
    if (method_bit > 0) assert(max_bits < method_bit);

    using namespace data_type;
    writer_b8_l ret(2);
    size_t adr = 0; size_t ofs_pos = 0;
    for (size_t curr = 0; adr < input.size(); ) {
      const auto& cmd = dp[curr][adr];
      const auto [tag, next, ulen] = cmd.type;
      switch (tag) {
      case uncomp: {
        ret.write<b1, d8>(false, input[adr]);
      } break;
      case lz: {
        const size_t lz_len = cmd.len - ulen;
        ret.write<b1, d16>(true, (adr - cmd.lz_ofs()) | (lz_len - lz_min_len) << (16 - len_bits));
      } break;
      default: assert(0);
      }
      if (next == 0 && (curr != 1 || ulen > 0)) {
        const size_t bits = 8 - ret.bit; assert(bits >= 1);
        const size_t bits_pos = ret.bits_pos;
        assert(ofs_pos + 2 <= ret.size());
        write16(ret.out, ofs_pos, bits_pos);
        ret.write<bnh>({ret.bit, ulen == 0 ? low_bits_mask(ret.bit) : 0});
        ret.write<d24>(0);
        for (size_t i = ret.size() - 1; i - 3 >= bits_pos; --i) ret[i] = ret[i - 3];
        if (ulen > 0) ret.write<d8n>({ulen, &input[adr + (cmd.len - ulen)]});
        ret[bits_pos + 0] = (bits + ulen) | method_bit;
        ofs_pos = bits_pos + 1;
      }
      adr += cmd.len;
      curr = next;
    }
    write16(ret.out, ofs_pos, ret.size());
    write16(ret.out, 0, read16(ret.out, 0) - 2);
    ret.write<d8>(method_bit);
    assert(adr == input.size());
    assert(dp[0].optimal_cost() + 3 == ret.size());

    if (best.empty() || ret.size() < best.size()) best = std::move(ret.out);
  }
  return best;
}


} // namespace

std::vector<uint8_t> bahamut_lagoon_comp(std::span<const uint8_t> input) {
  return chrono_trigger_comp_core(input, 4, 0x80);
}

std::vector<uint8_t> chrono_trigger_comp(std::span<const uint8_t> input) {
  return chrono_trigger_comp_core(input, 5, 0x3f);
}

std::vector<uint8_t> bahamut_lagoon_comp_fast(std::span<const uint8_t> input) {
  return chrono_trigger_comp_fast_core(input, 4);
}

std::vector<uint8_t> chrono_trigger_comp_fast(std::span<const uint8_t> input) {
  return chrono_trigger_comp_fast_core(input, 5);
}

} // namespace sfc_comp
