#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

struct rnc1_cost {
  using value_type = size_t;
  explicit constexpr rnc1_cost() : rnc1_cost(0, 0) {}
  explicit constexpr rnc1_cost(value_type value, size_t count = 0) : value(value), count(count) {}
  constexpr bool operator == (const rnc1_cost& rhs) const { return value == rhs.value && count == rhs.count; }
  constexpr bool operator != (const rnc1_cost& rhs) const { return !(*this == rhs); }
  constexpr rnc1_cost operator + (const value_type& c) const { return rnc1_cost(value + c, count); }
  constexpr bool operator >= (const rnc1_cost& rhs) const { return value >= rhs.value; }
  constexpr bool operator > (const rnc1_cost& rhs) const { return value > rhs.value; }
  constexpr bool operator < (const rnc1_cost& rhs) const { return !(*this >= rhs); }
  value_type value;
  size_t count; // uncomp count
};

} // namespace

template <>
struct cost_traits<rnc1_cost> {
  static constexpr rnc1_cost infinity() { return rnc1_cost(cost_traits<size_t>::infinity()); }
  static constexpr rnc1_cost unspecified() { return rnc1_cost(cost_traits<size_t>::unspecified()); }
};

std::vector<uint8_t> rob_northen_comp_1(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x100000);

  enum Tag { uncomp, lz };
  struct CompType {
    bool operator == (const CompType& rhs) const {
      if (tag != rhs.tag) return false;
      if (tag == uncomp) return li == rhs.li;
      return li == rhs.li;
    }
    Tag tag;
    size_t oi, li;
  };

  static constexpr size_t lz_min_len = 2;
  static constexpr size_t lz_ofs_max_bits = 15;
  static constexpr size_t lz_len_max_bits = 15; // <= 15
  static constexpr size_t uncomp_len_max_bits = 15; // <= 15

  static constexpr size_t command_limit = 0x4000;
  static_assert(1 <= command_limit && command_limit <= 0xffff, "invalid command_limit");

  static constexpr auto lz_offsets = create_array<vrange, lz_ofs_max_bits + 1>([](size_t k) {
    if (k == 0) return vrange(1, 1, 0, 0);
    return vrange((1 << (k - 1)) + 1, 1 << k, k - 1, 0);
  });
  static constexpr auto lz_lens = create_array<vrange, lz_len_max_bits + 1>([](size_t k) {
    if (k == 0) return vrange(lz_min_len, lz_min_len, 0, 0);
    return vrange((1 << (k - 1)) + lz_min_len, (1 << k) + (lz_min_len - 1), k - 1, 0);
  });
  static constexpr auto lz_max_len = lz_lens.back().max;
  static constexpr auto ulens = create_array<vrange, uncomp_len_max_bits + 1>([](size_t k) {
    if (k == 0) return vrange(0, 0, 0, 0);
    return vrange(1 << (k - 1), (1 << k) - 1, k - 1, 0);
  });

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, lz_ofs_max_bits + 1>> lz_memo(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    size_t o = lz_ofs_max_bits;
    auto res_lz = lz_helper.find_best_closest(i, lz_offsets[o].max, lz_max_len);
    if (res_lz.len < lz_min_len) res_lz = {0, 0};
    lz_memo[i][o] = res_lz;
    for (; o-- > 0; ) {
      const size_t d = i - res_lz.ofs;
      if (res_lz.len >= lz_min_len && d > lz_offsets[o].max) {
        res_lz = lz_helper.find_best_closest(i, lz_offsets[o].max, lz_max_len);
      }
      lz_memo[i][o] = res_lz;
    }
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16_l ret(0x12);
  ret.write<bnl>({2, 0});

  size_t chunk = 0;

  sssp_solver<CompType, rnc1_cost> dp0(input.size()), dp1(input.size());
  uncomp_helper<rnc1_cost::value_type> u_helper(input.size(), 8);

  size_t begin = 0;
  size_t last_index = begin;

  while (true) {
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> huff_u_bits(ulens.size());
    std::iota(huff_u_bits.begin(), huff_u_bits.end(), 1);

    std::vector<size_t> huff_lz_len_bits(lz_lens.size());
    std::iota(huff_lz_len_bits.begin(), huff_lz_len_bits.end(), 1);

    std::vector<size_t> huff_lz_ofs_bits(lz_offsets.size(), 4);

    struct huffman {
      encode::huffman_result uncomp_len;
      encode::huffman_result lz_len;
      encode::huffman_result lz_ofs;
    };

    auto update_costs = [&] (const encode::huffman_result& huff, std::vector<size_t>& costs, size_t penalty = 9) {
      size_t longest_bitlen = 0;
      for (const auto w : huff.words) longest_bitlen = std::max<size_t>(longest_bitlen, huff.codewords[w].bitlen);
      costs.assign(costs.size(), longest_bitlen + penalty);
      for (const auto w : huff.words) costs[w] = huff.codewords[w].bitlen;
    };

    auto update_huffman_costs = [&] (const huffman& huff) {
      update_costs(huff.uncomp_len, huff_u_bits);
      update_costs(huff.lz_ofs, huff_lz_ofs_bits);
      update_costs(huff.lz_len, huff_lz_len_bits);
    };

    using command_type = decltype(dp0)::vertex_type;
    size_t best_cost = std::numeric_limits<size_t>::max();
    size_t best_end = begin;
    std::vector<command_type> best_commands;
    huffman best_huff;

    while (true) {
      dp0[begin].cost = rnc1_cost(0);
      dp1.reset(begin);

      size_t e = std::min(input.size(), last_index + std::max(lz_lens.back().max, ulens.back().max));
      dp0.reset(begin + 1, e + 1);
      dp1.reset(begin + 1, e + 1);

      size_t index = begin;
      for (; index <= input.size(); ++index) {
        const auto cost0 = dp0[index].cost;
        if (cost0 < dp0.infinite_cost) {
          // skip uncomp
          const auto ncost = rnc1_cost(cost0.value + huff_u_bits[0], cost0.count + 1);
          dp1.update(index, 0, 0, Constant<0>(), {uncomp, 0, 0}, ncost);
        }
        if (index == input.size()) break;

        if (cost0 < dp0.infinite_cost) {
          u_helper.update(index, cost0.value);
        }

        for (size_t k = 1; k < ulens.size(); ++k) {
          const auto u = u_helper.find(index + 1, ulens[k].min, ulens[k].max);
          if (u.len == u_helper.nlen) continue;
          const auto cost = rnc1_cost(u.cost + huff_u_bits[k] + ulens[k].bitlen,
                                      dp0[index + 1 - u.len].cost.count + 1);
          dp1.update_u(index + 1, u.len, {uncomp, 0, k}, cost);
        }

        const auto cost1 = dp1[index].cost;
        if (cost1 == dp1.infinite_cost) break;
        if (cost1.count >= command_limit) break;

        for (ptrdiff_t oi = lz_ofs_max_bits; oi >= 0; --oi) {
          const auto& res_lz = lz_memo[index][oi];
          if (res_lz.len < lz_min_len) break;
          if ((index - res_lz.ofs) < lz_offsets[oi].min) continue;
          for (size_t li = 0; li < lz_lens.size(); ++li) {
            const size_t bit_size = (huff_lz_ofs_bits[oi] + lz_offsets[oi].bitlen) +
                                    (huff_lz_len_bits[li] + lz_lens[li].bitlen);
            const auto cost = rnc1_cost(cost1.value + bit_size, cost1.count);
            dp0.update_lz(index, lz_lens[li].min, lz_lens[li].max, res_lz,
                          Constant<0>(), {lz, size_t(oi), li}, cost);
          }
        }
      }
      u_helper.reset(begin, index);

      last_index = index;
      if (dp1[index].cost == dp1.infinite_cost) {
        if (index == 0) throw std::logic_error("This should not happen.");
        --index;
      }

      struct Counter {
        std::array<size_t, uncomp_len_max_bits + 1> uncomp = {};
        std::array<size_t, lz_len_max_bits + 1> lz_len = {};
        std::array<size_t, lz_ofs_max_bits + 1> lz_ofs = {};
      } counter;

      const auto commands = [&] {
        std::vector<command_type> ret;
        size_t adr = index;
        size_t ty = 1;
        while (adr > begin || (adr == begin && ty == 1)) {
          command_type cmd;
          if (ty == 1) cmd = dp1[adr], ty = 0;
          else cmd = dp0[adr], ty = 1;
          adr -= cmd.len;
          if (cmd.type.tag == uncomp) {
            counter.uncomp[cmd.type.li] += 1;
          } else {
            counter.lz_ofs[cmd.type.oi] += 1;
            counter.lz_len[cmd.type.li] += 1;
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
          ret.codewords[ret.words[0]].bitlen += 1;
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
        for (size_t w : huff.uncomp_len.words) {
          cost += counter.uncomp[w] * (huff.uncomp_len.codewords[w].bitlen + ulens[w].bitlen);
        }
        for (size_t w : huff.lz_ofs.words) {
          cost += counter.lz_ofs[w] * (huff.lz_ofs.codewords[w].bitlen + lz_offsets[w].bitlen);
        }
        for (size_t w : huff.lz_len.words) {
          cost += counter.lz_len[w] * (huff.lz_len.codewords[w].bitlen + lz_lens[w].bitlen);
        }
        return cost;
      };

      const size_t estimated_cost = calc_cost(huff, commands, counter);

      if (estimated_cost < best_cost) {
        best_cost = estimated_cost;
        best_commands = std::move(commands);
        best_huff = std::move(huff);
        best_end = index; // (do not use last_index)
        update_huffman_costs(best_huff);
      } else {
        auto write_cost = [&](const encode::huffman_result& huff) {
          const size_t total_count = max_elem(huff.words) + 1;
          ret.write<bnl>({5, total_count});
          for (size_t i = 0; i < total_count; ++i) {
            ptrdiff_t bits = huff.codewords[i].bitlen;
            if (bits <= 0) bits = 0;
            if (bits >= 16) {
              // this would not happen.
              throw std::runtime_error("Failed to compress the given data.");
            }
            ret.write<bnl>({4, size_t(bits)});
          }
        };
        write_cost(best_huff.uncomp_len);
        write_cost(best_huff.lz_ofs);
        write_cost(best_huff.lz_len);
        assert(best_commands.size() % 2 == 1 && (best_commands.size() + 1) / 2 <= command_limit);
        ret.write<bnl>({16, (best_commands.size() + 1) / 2});

        size_t adr = begin;
        for (const auto& cmd : best_commands) {
          switch (cmd.type.tag) {
          case uncomp: {
            const size_t k = cmd.type.li;
            ret.write<bnh>(best_huff.uncomp_len.codewords[k]);
            if (cmd.len > 0) {
              assert(ulens[k].min <= cmd.len && cmd.len <= ulens[k].max);
              ret.write<bnl>({ulens[k].bitlen, cmd.len - ulens[k].min});
              ret.write<d8n>({cmd.len, &input[adr]});
            }
          } break;
          case lz: {
            const size_t oi = cmd.type.oi;
            const auto& o = lz_offsets[oi];
            const size_t d = (adr - cmd.lz_ofs);
            assert(o.min <= d && d <= o.max);
            ret.write<bnh>(best_huff.lz_ofs.codewords[oi]);
            ret.write<bnl>({o.bitlen, d - o.min});

            const size_t li = cmd.type.li;
            const auto l = lz_lens[li];
            assert(l.min <= cmd.len && cmd.len <= l.max);
            ret.write<bnh>(best_huff.lz_len.codewords[li]);
            ret.write<bnl>({l.bitlen, cmd.len - l.min});
          } break;
          default: assert(0);
          }
          adr += cmd.len;
        }
        assert(adr == best_end);
        begin = best_end;
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
