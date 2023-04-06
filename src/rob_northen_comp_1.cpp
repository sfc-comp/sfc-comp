#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

struct CostType {
  explicit constexpr CostType() : CostType(0, 0) {}
  explicit constexpr CostType(size_t cost, size_t count = 0) : cost(cost), count(count) {}

  constexpr bool operator == (const CostType& rhs) const { return cost == rhs.cost && count == rhs.count; }
  constexpr bool operator != (const CostType& rhs) const { return !(*this == rhs); }
  constexpr CostType operator + (const CostType& rhs) const {
    return CostType(cost + rhs.cost, count + rhs.count);
  }
  constexpr CostType operator + (const size_t& c) const {
    return CostType(cost + c, count);
  }
  constexpr bool operator >= (const CostType& rhs) const {
    return cost >= rhs.cost;
  }
  constexpr bool operator < (const CostType& rhs) const {
    return !(*this >= rhs);
  }
  size_t cost;
  size_t count; // uncomp count
};

} // namespace

template <>
constexpr CostType default_cost() { return CostType(default_cost<size_t>()); }

std::vector<uint8_t> rob_northen_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x100000);

  enum Tag {
    uncomp, lz
  };

  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return len_no == rhs.len_no;
      return len_no == rhs.len_no && ofs_no == rhs.ofs_no;
    }
    Tag tag;
    size_t ofs_no, len_no;
  };

  static constexpr size_t lz_max_ofs_bits = 15;
  static constexpr size_t lz_max_len_bits = 15; // <= 15
  static constexpr size_t uncomp_max_len_bits = 10; // <= 15

  static constexpr size_t command_limit = 0x4000;
  static_assert(1 <= command_limit && command_limit <= 0xffff, "invalid command_limit");

  std::array<size_t, lz_max_ofs_bits + 1> lz_ofs_table = {};
  for (size_t i = 0; i < lz_ofs_table.size(); ++i) lz_ofs_table[i] = 1 << i;
  std::array<size_t, lz_max_len_bits + 1> lz_len_table = {};
  for (size_t i = 0; i < lz_len_table.size(); ++i) lz_len_table[i] = (1 << i) + 1;
  std::array<size_t, uncomp_max_len_bits + 1> uncomp_len_table = {};
  for (size_t i = 0; i < uncomp_len_table.size(); ++i) uncomp_len_table[i] = (1 << i) - 1;

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, lz_max_ofs_bits + 1>> lz_memo(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    auto res_lz = lz_helper.find_best_closest(i, lz_ofs_table.back(), lz_len_table.back());
    if (res_lz.len < lz_len_table[0]) res_lz = {0, 0};
    size_t o = lz_max_ofs_bits;
    lz_memo[i][o] = res_lz;
    for (; o-- > 0; ) {
      const size_t d = i - res_lz.ofs;
      if (res_lz.len >= lz_len_table[0] && d > lz_ofs_table[o]) {
        res_lz = lz_helper.find_best_closest(i, lz_ofs_table[o], lz_len_table.back());
      }
      lz_memo[i][o] = res_lz;
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16 ret;
  for (size_t i = 0; i < 0x12; ++i) ret.write<d8>(0);
  ret.write<b8ln_l>({2, 0});

  size_t chunk = 0;

  using command_type = sssp_solver<CompType, CostType>::vertex_type;
  sssp_solver<CompType, CostType> dp0(input.size());
  sssp_solver<CompType, CostType> dp1(input.size());

  size_t begin = 0;
  size_t last_index = begin;

  while (true) {
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> huff_uncomp_bits(uncomp_len_table.size());
    std::iota(huff_uncomp_bits.begin(), huff_uncomp_bits.end(), 1);

    std::vector<size_t> huff_lz_len_bits(lz_len_table.size());
    std::iota(huff_lz_len_bits.begin(), huff_lz_len_bits.end(), 1);

    std::vector<size_t> huff_lz_ofs_bits(lz_ofs_table.size());
    std::iota(huff_lz_ofs_bits.begin(), huff_lz_ofs_bits.end(), 1);

    struct huffman {
      encode::huffman_result uncomp_len;
      encode::huffman_result lz_len;
      encode::huffman_result lz_ofs;
    };

    auto update_costs = [&] (const encode::huffman_result& huff, std::vector<size_t>& costs, size_t penalty = 9) {
      size_t worst_bit_count = 0;
      for (const auto w : huff.words) worst_bit_count = std::max<size_t>(worst_bit_count, huff.codewords[w].bit_count);
      costs.assign(costs.size(), worst_bit_count + penalty);
      for (const auto w : huff.words) costs[w] = huff.codewords[w].bit_count;
    };

    auto update_huffman_costs = [&] (const huffman& huff) {
      update_costs(huff.uncomp_len, huff_uncomp_bits);
      update_costs(huff.lz_ofs, huff_lz_ofs_bits);
      update_costs(huff.lz_len, huff_lz_len_bits);
    };

    size_t best_cost = std::numeric_limits<size_t>::max();
    size_t best_last_index = 0;
    std::vector<command_type> best_commands;
    huffman best_huff;

    while (true) {
      dp0[begin].cost = CostType(0);
      dp1[begin].cost = dp0.default_cost;

      size_t e = std::min(input.size(), last_index + std::max(lz_len_table.back(), uncomp_len_table.back()));
      for (size_t i = begin + 1; i <= e; ++i) {
        dp0[i].cost = dp0.default_cost;
        dp1[i].cost = dp1.default_cost;
      }

      constexpr auto bit_cost = [](size_t b) -> size_t {
        return (b == 0) ? 0 : (b - 1);
      };

      size_t index = begin;
      for (; index <= input.size(); ++index) {
        const auto cost0 = dp0[index].cost;
        if (cost0 != dp0.default_cost) {
          // skip uncomp
          const auto ncost = CostType(cost0.cost + huff_uncomp_bits[0], cost0.count + 1);
          if (ncost < dp1[index].cost) {
            dp1[index].cost = ncost;
            dp1[index].len = 0;
            dp1[index].lz_ofs = 0;
            dp1[index].type = {uncomp, 0, 0};
          }
        }
        if (index == input.size()) break;

        if (cost0 != dp0.default_cost) {
          // uncomp
          for (size_t k = 1; k < uncomp_len_table.size(); ++k) {
            const size_t from = uncomp_len_table[k - 1] + 1;
            const size_t to = uncomp_len_table[k];
            const auto base_cost = CostType(cost0.cost + huff_uncomp_bits[k] + bit_cost(k), cost0.count + 1);
            dp1.update(index, from, to, Linear<8, 0>(), {uncomp, 0, k}, base_cost);
          }
        }

        const auto cost1 = dp1[index].cost;
        if (cost1 == dp1.default_cost) break;
        if (cost1.count >= command_limit) break;

        for (ptrdiff_t o = lz_max_ofs_bits; o >= 0; --o) {
          const auto& res_lz = lz_memo[index][o];
          if (res_lz.len < lz_len_table[0]) break;
          const size_t lz_min_dist = (o == 0) ? lz_ofs_table[0] : (lz_ofs_table[o - 1] + 1);
          if ((index - res_lz.ofs) < lz_min_dist) continue;
          for (size_t k = 0; k < lz_len_table.size(); ++k) {
            const size_t from = (k == 0) ? lz_len_table[0] : (lz_len_table[k - 1] + 1);
            const size_t to = lz_len_table[k];
            const size_t bit_size = (huff_lz_ofs_bits[o] + bit_cost(o)) +
                                    (huff_lz_len_bits[k] + bit_cost(k));
            const auto base_cost = CostType(cost1.cost + bit_size, cost1.count);
            dp0.update_lz(index, from, to, res_lz, Constant<0>(), {lz, size_t(o), k}, base_cost);
          }
        }
      }

      last_index = index;
      if (dp1[index].cost == dp1.default_cost) {
        if (index == 0) throw std::logic_error("This should not happen.");
        --index;
      }

      struct Counter {
        Counter() : uncomp(uncomp_max_len_bits + 1, 0),
                    lz_len(lz_max_len_bits + 1, 0),
                    lz_ofs(lz_max_ofs_bits + 1, 0) {}
        std::vector<size_t> uncomp;
        std::vector<size_t> lz_len;
        std::vector<size_t> lz_ofs;
      } counter;

      const std::vector<command_type> commands = [&] {
        std::vector<command_type> ret;
        size_t adr = index;
        size_t ty = 1;
        while (adr > begin || (adr == begin && ty == 1)) {
          command_type cmd;
          if (ty == 1) cmd = dp1[adr], ty = 0;
          else cmd = dp0[adr], ty = 1;
          adr -= cmd.len;
          if (cmd.type.tag == uncomp) {
            counter.uncomp[cmd.type.len_no] += 1;
          } else {
            counter.lz_ofs[cmd.type.ofs_no] += 1;
            counter.lz_len[cmd.type.len_no] += 1;
          }
          ret.emplace_back(cmd);
        }
        assert(adr == begin && ty == 0);
        std::reverse(ret.begin(), ret.end());
        return ret;
      }();

      auto huffman_encode = [&](std::span<const size_t> counts) {
        auto ret = encode::huffman(counts, true);
        if (ret.words.size() == 1) {
          // 0-bit codes are not allowed.
          ret.codewords[ret.words[0]].bit_count += 1;
        }
        return ret;
      };

      huffman huff;
      huff.uncomp_len = huffman_encode(counter.uncomp);
      huff.lz_ofs = huffman_encode(counter.lz_ofs);
      huff.lz_len = huffman_encode(counter.lz_len);

      auto max_elem = [&](std::span<const size_t> words) -> ptrdiff_t {
        if (words.size() == 0) return -1; // This value must be -1.
        return *std::max_element(words.begin(), words.end());
      };

      auto calc_cost = [&](const huffman& huff, std::span<const command_type> commands, const Counter& counter) -> size_t {
        size_t cost = 0;
        cost += 5 + 4 * (1 + max_elem(huff.uncomp_len.words));
        cost += 5 + 4 * (1 + max_elem(huff.lz_ofs.words));
        cost += 5 + 4 * (1 + max_elem(huff.lz_len.words));
        cost += 16;
        for (const auto& cmd : commands) {
          if (cmd.type.tag == uncomp) cost += cmd.len * 8;
        }
        for (size_t w : huff.uncomp_len.words) cost += counter.uncomp[w] * (huff.uncomp_len.codewords[w].bit_count + bit_cost(w));
        for (size_t w : huff.lz_ofs.words) cost += counter.lz_ofs[w] * (huff.lz_ofs.codewords[w].bit_count + bit_cost(w));
        for (size_t w : huff.lz_len.words) cost += counter.lz_len[w] * (huff.lz_len.codewords[w].bit_count + bit_cost(w));
        return cost;
      };

      const size_t estimated_cost = calc_cost(huff, commands, counter);

      if (estimated_cost < best_cost) {
        best_cost = estimated_cost;
        best_commands = std::move(commands);
        best_huff = std::move(huff);
        best_last_index = index; // (do not use last_index)
        update_huffman_costs(best_huff);
      } else {
        auto write_cost = [&](const encode::huffman_result& huff) {
          const size_t total_count = max_elem(huff.words) + 1;
          ret.write<b8ln_l>({5, total_count});
          for (size_t i = 0; i < total_count; ++i) {
            ptrdiff_t bits = huff.codewords[i].bit_count;
            if (bits <= 0) bits = 0;
            if (bits >= 16) {
              // this would not happen.
              throw std::runtime_error("Failed to compress the given data.");
            }
            ret.write<b8ln_l>({4, size_t(bits)});
          }
        };
        write_cost(best_huff.uncomp_len);
        write_cost(best_huff.lz_ofs);
        write_cost(best_huff.lz_len);
        assert(best_commands.size() % 2 == 1);
        ret.write<b8ln_l>({16, (best_commands.size() + 1) / 2});

        size_t adr = begin;
        for (const auto cmd : best_commands) {
          switch (cmd.type.tag) {
          case uncomp: {
            const size_t k = cmd.type.len_no;
            const auto c = best_huff.uncomp_len.codewords[k];
            assert(c.bit_count >= 0);
            ret.write<b8ln_h>({size_t(c.bit_count), c.val});
            if (cmd.len > 0) {
              const size_t min_len = ((k >= 1) ? uncomp_len_table[k - 1] + 1 : uncomp_len_table[0]);
              assert(min_len <= cmd.len && cmd.len <= uncomp_len_table[k]);
              if (k >= 2) ret.write<b8ln_l>({k - 1, cmd.len - min_len});
              ret.write<d8n>({cmd.len, &input[adr]});
            }
          } break;
          case lz: {
            const size_t k = cmd.type.len_no, o = cmd.type.ofs_no;
            const auto c_len = best_huff.lz_len.codewords[k];
            const auto c_ofs = best_huff.lz_ofs.codewords[o];

            const size_t d = (adr - cmd.lz_ofs);
            const size_t min_ofs = ((o >= 1) ? lz_ofs_table[o - 1] + 1 : lz_ofs_table[0]);
            assert(min_ofs <= d && d <= lz_ofs_table[o]);
            ret.write<b8ln_h>({size_t(c_ofs.bit_count), c_ofs.val});
            if (o >= 2) ret.write<b8ln_l>({o - 1, d - min_ofs});

            const size_t min_len = ((k >= 1) ? lz_len_table[k - 1] + 1 : lz_len_table[0]);
            assert(min_len <= cmd.len && cmd.len <= lz_len_table[k]);
            ret.write<b8ln_h>({size_t(c_len.bit_count), c_len.val});
            if (k >= 2) ret.write<b8ln_l>({k - 1, cmd.len - min_len});
          } break;
          default: assert(0);
          }
          adr += cmd.len;
        }
        assert(adr == best_last_index);
        begin = adr;
        break;
      }
    }
    ++chunk;
    if (chunk == 0x100) throw std::runtime_error("This algorithm cannot compress the given data. (chunk >= 0x100)");
    if (begin == input.size()) break;
  }

  write32b(ret.out, 0, 0x524e4301);
  write32b(ret.out, 4, input.size());
  write32b(ret.out, 8, ret.size() - 18);
  write16b(ret.out, 12, utility::crc16(input, 0, input.size()));
  write16b(ret.out, 14, utility::crc16(ret.out, 18, ret.size() - 18));
  ret.out[0x10] = 0; // leeway
  ret.out[0x11] = chunk;

  return ret.out;
}

} // namespace sfc_comp
