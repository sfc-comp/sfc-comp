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

    lz_helper lz_helper(input);
    sssp_solver<tag> dp(input.size());

    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, Constant<9>(), uncomp);
      auto res_lz = lz_helper.find_best(i, lz_max_ofs);
      dp.update_lz(i, 3, lz_max_len, res_lz, Constant<17>(), lz);
      lz_helper.add_element(i);
    }

    using namespace data_type;
    writer_b8_l ret(2);
    size_t adr = 0;
    for (const auto cmd : dp.commands()) {
      switch (cmd.type) {
      case uncomp: ret.write<b1, d8>(false, input[adr]); break;
      case lz: ret.write<b1, d16>(true, (adr - cmd.lz_ofs) | (cmd.len - lz_min_len) << (16 - len_bits)); break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    assert(adr == input.size());
    assert(dp.total_cost() + 2 * 8 == ret.bit_length());

    const size_t method_bit = (len_bits == 4) ? 0x00 : 0x40;
    if (ret.bit == 0) {
      write16(ret.out, 0, ret.size() - 2);
    } else {
      const size_t bits_pos = ret.bits_pos;
      write16(ret.out, 0, bits_pos - 2);

      const size_t len = (8 - ret.bit) + popcount32(ret.out[bits_pos]) + 1;
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

  enum method { uncomp, lz, none };
  using tag = tag_ol<method>;

  std::vector<uint8_t> best;

  for (size_t len_bits = 4; len_bits <= max_len_bits; ++len_bits) {
    const size_t lz_min_len = 3;
    const size_t lz_max_len = ((1 << len_bits) - 1) + lz_min_len;
    const size_t lz_max_ofs = (0x10000 >> len_bits) - 1;

    lz_helper lz_helper(input);
    std::array<sssp_solver<tag>, 8> dp;
    for (size_t b = 0; b < 8; ++b) dp[b] = sssp_solver<tag>(input.size(), b == 0 ? 0 : -1);

    for (size_t i = 0; i < input.size(); ++i) {
      const auto res_lz = lz_helper.find_best(i, lz_max_ofs);
      const auto update = [&](size_t b, size_t from, auto cost) {
        dp[b].update(i, 1, 1, Constant<1>(), {uncomp, from, 0}, cost);
        dp[b].update_lz(i, lz_min_len, lz_max_len, res_lz, Constant<2>(), {lz, from, 0}, cost);
      };
      for (size_t b = 8; b-- > 0; ) {
        const auto cost = dp[b][i].cost;
        if (cost == dp[0].infinite_cost) continue;
        const auto ncost = cost + (b == 0 ? 1 : 0);
        update((b - 1) & 7, b, ncost);
        if (b != 1) update(0, b, ncost + 3);
        if (res_lz.len >= lz_min_len) {
          const size_t remain = max_bits - 1 - ((8 - b) & 7);
          const size_t lz_len = std::min(res_lz.len, lz_max_len);
          for (size_t u = 1; u <= remain; ++u) {
            dp[0].update(i, lz_len + u, lz_len + u, Constant<2>(), {lz, b, u}, ncost + 3 + u, res_lz.ofs);
          }
        }
      }
      lz_helper.add_element(i);
    }

    const auto [commands, offset] = [&] {
      using command_type = sssp_solver<tag>::vertex_type;
      std::vector<command_type> ret;
      ptrdiff_t adr = input.size();
      size_t curr = 0; size_t offset = dp[curr][adr].cost, bits = 0;
      while (adr > 0 || (adr == 0 && curr > 0)) {
        auto cmd = dp[curr][adr];
        const size_t ulen = cmd.type.li;
        const size_t prev = cmd.type.oi;
        if (curr == 0 && (prev != 1 || ulen > 0)) {
          bits = ((curr - prev) & 7) + 1 + ulen;
        }
        curr = prev;
        adr -= cmd.len;
        ret.emplace_back(cmd);

        if (curr == 0 && bits > 0) {
          ret.emplace_back(-1, 0, offset, tag(none, 0, bits));
          offset = dp[curr][adr].cost; bits = 0;
        }
      }
      assert(adr == 0 && curr == 0);
      std::reverse(ret.begin(), ret.end());
      return std::make_pair(std::move(ret), offset);
    }();

    const size_t method_bit = (len_bits == 4) ? 0x00 : 0x40;
    if (method_bit > 0) assert(max_bits < method_bit);

    using namespace data_type;
    writer_b8_l ret(2);

    size_t bits = 0;
    const auto fix = [&](size_t b, bool avoid_zero) {
      if (bits == 0) return;
      assert(bits >= b); bits -= b;
      if (bits == 0) {
        if (avoid_zero) ret.write<bnh>({ret.bit, low_bits_mask(ret.bit)});
        ret.bit = 0;
      }
    };

    size_t adr = 0;
    for (const auto& cmd : commands) {
      switch (cmd.type.tag) {
      case none: {
        assert(ret.bit == 0);
        bits = cmd.type.li; assert(bits <= max_bits);
        ret.write<d8, d16>(bits | method_bit, cmd.val() + 2);
      } break;
      case uncomp: {
        ret.write<b1, d8>(false, input[adr]); fix(1, true);
      } break;
      case lz: {
        const size_t ulen = cmd.type.li;
        const size_t lz_len = cmd.len - ulen;
        ret.write<b1, d16>(true, (adr - cmd.lz_ofs) | (lz_len - lz_min_len) << (16 - len_bits)); fix(1, true);
        if (ulen > 0) {
          ret.write<d8n>({ulen, &input[adr + lz_len]}); fix(ulen, false);
        }
      } break;
      default: assert(0);
      }
      adr += cmd.len;
    }
    write16(ret.out, 0, offset);
    ret.write<d8>(method_bit);
    assert(adr == input.size());
    assert(dp[0].total_cost() + 3 == ret.size());

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
