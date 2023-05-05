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
      ret[bits_pos + 3] |= 1 << (8 - ret.bit); // Avoids 0x00. (cf. $C3:07B7, $C3:0879, etc. in Chrono Trigger)
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
    std::array<sssp_solver<tag>, 16> dp;
    for (size_t b = 0; b < 16; ++b) dp[b] = sssp_solver<tag>(input.size(), b == 0 ? 0 : -1);

    for (size_t i = 0; i < input.size(); ++i) {
      const auto res_lz = lz_helper.find_best(i, lz_max_ofs);
      const auto update = [&](size_t b, size_t from, auto cost) {
        dp[b].update(i, 1, 1, Constant<1>(), {uncomp, from, 0}, cost);
        dp[b].update_lz(i, lz_min_len, lz_max_len, res_lz, Constant<2>(), {lz, from, 0}, cost);
      };
      for (size_t b = 7; b > 0; --b) {
        const auto cost = dp[b][i].cost;
        if (cost == dp[0].infinite_cost) continue;
        update(b - 1, b, cost);
      }
      if (const auto cost = dp[0][i].cost; cost < dp[0].infinite_cost) {
        update(7, 0, cost + 1);
        dp[15].update(i, 0, 0, Constant<1 + 2>(), {none, 0, 0}, cost + 1);
      }
      for (size_t b = 15; b >= 8; --b) {
        const auto cost = dp[b][i].cost;
        if (cost == dp[0].infinite_cost) continue;
        update(0, b, cost);
        if (b > 8) update(b - 1, b, cost);
        const size_t remain = max_bits - 1 - (15 - b);
        if (res_lz.len >= lz_min_len) {
          // lz + uncomp
          const size_t lz_len = std::min(res_lz.len, lz_max_len);
          for (size_t u = 1; u <= remain; ++u) {
            dp[0].update(i, lz_len + u, lz_len + u, Constant<2>(), {lz, b, u}, cost + u, res_lz.ofs);
          }
        }
      }
      lz_helper.add_element(i);
    }

    const auto [commands, offset] = [&] {
      using command_type = sssp_solver<tag>::vertex_type;
      std::vector<command_type> ret;
      ptrdiff_t adr = input.size();
      size_t curr = 0; size_t offset = 0, bits = 0;
      while (adr > 0 || (adr == 0 && curr > 0)) {
        auto cmd = dp[curr][adr];
        curr = cmd.type.oi; // previous state
        if (curr >= 8) bits += 1 + cmd.type.li, offset += cmd.type.li;
        if (curr == 0) offset += 1;
        const auto tag = cmd.type.tag;
        if (tag == uncomp) {
          offset += 1;
        } else if (tag == lz) {
          offset += 2;
        } else {
          cmd.set_val(offset); offset = 0;
          cmd.type.li = bits; bits = 0;
        }
        adr -= cmd.len;
        ret.emplace_back(cmd);
      }
      assert(adr == 0 && curr == 0);
      std::reverse(ret.begin(), ret.end());
      return std::make_pair(std::move(ret), offset);
    }();

    const size_t method_bit = (len_bits == 4) ? 0x00 : 0x40;
    if (method_bit > 0) assert(max_bits < method_bit);

    using namespace data_type;
    writer_b8_l ret(2);

    size_t remaining_bits = 0;
    const auto fix = [&](size_t b, bool avoid_zero) {
      if (remaining_bits == 0) return;
      assert(remaining_bits >= b);
      remaining_bits -= b;
      if (remaining_bits == 0) {
        if (ret.bit > 0 && avoid_zero) ret.write<b1>(true);
        ret.bit = 0;
      }
    };

    size_t adr = 0;
    for (const auto& cmd : commands) {
      switch (cmd.type.tag) {
      case none: {
        const size_t s = ret.size();
        const size_t o = cmd.val();
        remaining_bits = cmd.type.li;
        assert(ret.bit == 0);
        assert(remaining_bits <= max_bits);
        ret.write<d8, d16>(remaining_bits | method_bit, s + 3 + o);
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
