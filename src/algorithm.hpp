#pragma once

#include <cstdio>
#include <cstdint>
#include <cstddef>
#include <ctime>
#include <cassert>

#include <vector>
#include <utility>
#include <limits>

#include <array>
#include <tuple>
#include <type_traits>

#include <span>

#include "encode.hpp"
#include "data_structure.hpp"

namespace sfc_comp {

namespace encode {

struct lz_data {
  lz_data() : ofs(0), len(0) {}
  lz_data(ptrdiff_t ofs, size_t len) : ofs(ofs), len(len) {}
  ptrdiff_t ofs;
  size_t len;
};

namespace lz {

template <typename SegtreeLcp, typename Segtree>
encode::lz_data find_best(
    size_t pos, size_t index, size_t block_size, const SegtreeLcp& seg_lcp, const Segtree& seg) {
  auto r = seg.find_range(index, [&](int v) {
    return v + static_cast<int>(block_size) < static_cast<int>(pos);
  });
  size_t len = 0; int offset = -1;
  if (r.first > 0) {
    size_t len_l = seg_lcp.fold(r.first - 1, index);
    if (len_l > len) len = len_l, offset = seg[r.first - 1];
  }
  if (r.second < seg.size()) {
    size_t len_r = seg_lcp.fold(index, r.second);
    if (len_r > len) len = len_r, offset = seg[r.second];
  }
  return {offset, len};
}

template <typename SegtreeLcp, typename Segtree>
encode::lz_data find_best_closest(
    size_t pos, size_t index, size_t block_size, size_t max_len,
    const SegtreeLcp& seg_lcp, const Segtree& seg) {
  auto ret = find_best(pos, index, block_size, seg_lcp, seg);
  if (ret.len > 0) {
    if (ret.len > max_len) ret.len = max_len;
    auto r2 = seg_lcp.find_range(index, [&](size_t v) { return v >= ret.len; });
    ret.ofs = seg.fold(r2.first, r2.second + 1);
  }
  return ret;
}

} // namespace lz

} // namespace encode

template <typename T>
struct range_max {
  using value_type = T;
  static T unit() { return std::numeric_limits<T>::min(); }
  static T op(const value_type& l, const value_type& r) { return std::max<T>(l, r); }
};

template <typename T>
struct range_min {
  using value_type = T;
  static T unit() { return std::numeric_limits<T>::max(); }
  static T op(const value_type& l, const value_type& r) { return std::min<T>(l, r); }
};

class lz_helper {
public:
  lz_helper(std::span<const uint8_t> input)
    : n(input.size()),
      segtree(n) {

    const auto plcp = suffix_array<uint8_t>(input).compute_lcp_and_rank(input);
    const auto& lcp = plcp.first;
    rank = std::move(plcp.second);
    segtree_lcp = segment_tree<range_min<int>>(lcp);
  }

  encode::lz_data find_best(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank[pos], block_size, segtree_lcp, segtree);
  }

  encode::lz_data find_best_closest(size_t pos, size_t block_size, size_t max_len) const {
    return encode::lz::find_best_closest(pos, rank[pos], block_size, max_len, segtree_lcp, segtree);
  }

  void add_element(size_t i) {
    segtree.update(rank[i], i);
  }

private:
  const size_t n;
  std::vector<int> rank;
  segment_tree<range_max<int>> segtree;
  segment_tree<range_min<int>> segtree_lcp;
};

class lz_helper_c {
private:
  std::vector<int16_t> complement_appended(std::span<const uint8_t> input) const {
    std::vector<int16_t> input_xor(2 * n + 1);
    for (size_t i = 0; i < n; ++i) input_xor[i] = input[i];
    input_xor[n] = -1;
    for (size_t i = 0; i < n; ++i) input_xor[i + n + 1] = input[i] ^ 0xff;
    return input_xor;
  }

public:
  lz_helper_c(std::span<const uint8_t> input)
      : n(input.size()) {

    const auto input_c = complement_appended(input);
    const auto plcp = suffix_array<int16_t>(input_c).compute_lcp_and_rank(input_c);
    const auto& lcp = plcp.first;
    rank = std::move(plcp.second);
    segtree = segment_tree<range_max<int>>(input_c.size());
    segtree_c = segment_tree<range_max<int>>(input_c.size());
    segtree_lcp = segment_tree<range_min<int>>(lcp);
  }

public:
  encode::lz_data find_best(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank[pos], block_size, segtree_lcp, segtree);
  }

  encode::lz_data find_best_c(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank[pos], block_size, segtree_lcp, segtree_c);
  }

  encode::lz_data find_best_closest(size_t pos, size_t block_size, size_t max_len) const {
    return encode::lz::find_best_closest(pos, rank[pos], block_size, max_len, segtree_lcp, segtree);
  }

  encode::lz_data find_best_closest_c(size_t pos, size_t block_size, size_t max_len) const {
    return encode::lz::find_best_closest(pos, rank[pos], block_size, max_len, segtree_lcp, segtree_c);
  }

  void add_element(size_t i) {
    segtree.update(rank[i], i);
    segtree_c.update(rank[i + n + 1], i);
  }

private:
  const size_t n;
  std::vector<int> rank;
  segment_tree<range_max<int>> segtree;
  segment_tree<range_max<int>> segtree_c;
  segment_tree<range_min<int>> segtree_lcp;
};

constexpr uint8_t bit_reversed[256] = {
  0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0,
  0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0, 0x70, 0xf0,
  0x08, 0x88, 0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8,
  0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8,
  0x04, 0x84, 0x44, 0xc4, 0x24, 0xa4, 0x64, 0xe4,
  0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4,
  0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac, 0x6c, 0xec,
  0x1c, 0x9c, 0x5c, 0xdc, 0x3c, 0xbc, 0x7c, 0xfc,
  0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2,
  0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2,
  0x0a, 0x8a, 0x4a, 0xca, 0x2a, 0xaa, 0x6a, 0xea,
  0x1a, 0x9a, 0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa,
  0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6,
  0x16, 0x96, 0x56, 0xd6, 0x36, 0xb6, 0x76, 0xf6,
  0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee,
  0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe, 0x7e, 0xfe,
  0x01, 0x81, 0x41, 0xc1, 0x21, 0xa1, 0x61, 0xe1,
  0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
  0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9,
  0x19, 0x99, 0x59, 0xd9, 0x39, 0xb9, 0x79, 0xf9,
  0x05, 0x85, 0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5,
  0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5,
  0x0d, 0x8d, 0x4d, 0xcd, 0x2d, 0xad, 0x6d, 0xed,
  0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd,
  0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3, 0x63, 0xe3,
  0x13, 0x93, 0x53, 0xd3, 0x33, 0xb3, 0x73, 0xf3,
  0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb,
  0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb,
  0x07, 0x87, 0x47, 0xc7, 0x27, 0xa7, 0x67, 0xe7,
  0x17, 0x97, 0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7,
  0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef,
  0x1f, 0x9f, 0x5f, 0xdf, 0x3f, 0xbf, 0x7f, 0xff
};

class lz_helper_kirby {
private:
  std::vector<int16_t> hflip_appended(std::span<const uint8_t> input) const {
    std::vector<int16_t> ret(2 * n + 1);
    for (size_t i = 0; i < n; ++i) ret[i] = input[i];
    ret[n] = -1;
    for (size_t i = 0; i < n; ++i) ret[i + n + 1] = bit_reversed[input[i]];
    return ret;
  }

  std::vector<int16_t> vflip_appended(std::span<const uint8_t> input) const {
    std::vector<int16_t> ret(2 * n + 1);
    for (size_t i = 0; i < n; ++i) ret[i] = input[i];
    ret[n] = -1;
    for (size_t i = 0; i < n; ++i) ret[i + n + 1] = input[n - 1 - i];
    return ret;
  }

public:
  lz_helper_kirby(std::span<const uint8_t> input)
      : n(input.size()) {

    const auto input_h = hflip_appended(input);
    const auto plcp_h = suffix_array<int16_t>(input_h).compute_lcp_and_rank(input_h);
    seg_lcp_h = segment_tree<range_min<int>>(plcp_h.first);
    rank_h = std::move(plcp_h.second);
    seg = segment_tree<range_max<int>>(input_h.size());
    seg_h = segment_tree<range_max<int>>(input_h.size());

    const auto input_v = vflip_appended(input);
    const auto plcp_v = suffix_array<int16_t>(input_v).compute_lcp_and_rank(input_v);
    seg_lcp_v = segment_tree<range_min<int>>(plcp_v.first);
    rank_v = std::move(plcp_v.second);
    seg_v = segment_tree<range_max<int>>(input_v.size());
  }

public:
  encode::lz_data find_best(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank_h[pos], block_size, seg_lcp_h, seg);
  }

  encode::lz_data find_best_h(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank_h[pos], block_size, seg_lcp_h, seg_h);
  }

  encode::lz_data find_best_v(size_t pos, size_t block_size) const {
    return encode::lz::find_best(pos, rank_v[pos], block_size, seg_lcp_v, seg_v);
  }

  void add_element(size_t i) {
    seg.update(rank_h[i], i);
    seg_h.update(rank_h[i + n + 1], i);
    seg_v.update(rank_v[2 * n - i], i);
  }

private:
  const size_t n;
  std::vector<int> rank_h;
  std::vector<int> rank_v;
  segment_tree<range_min<int>> seg_lcp_h;
  segment_tree<range_min<int>> seg_lcp_v;
  segment_tree<range_max<int>> seg;
  segment_tree<range_max<int>> seg_h;
  segment_tree<range_max<int>> seg_v;
};

template <size_t C>
struct Constant {
  size_t operator() (size_t) const {
    return C;
  }
};

template <size_t A, size_t B>
struct Linear {
  size_t operator() (size_t i) const {
    return A * i + B;
  }
};

template <size_t A, size_t B, size_t C>
struct LinearQ {
  size_t operator() (size_t i) const {
    return (A * i + B) / C;
  }
};

template <typename>
struct is_constant : std::false_type {};

template <size_t N>
struct is_constant<Constant<N>> : std::true_type {};

template <typename>
struct is_linear : std::false_type {};

template <size_t A, size_t B>
struct is_linear<Linear<A, B>> : std::true_type {};

template <typename, size_t>
struct is_linear_k : std::false_type {};

template <size_t A, size_t B, size_t K>
struct is_linear_k<Linear<A, B>, K> : std::true_type {};

template <size_t A, size_t B, size_t C, size_t K>
struct is_linear_k<LinearQ<A, B, C>, K>
  : std::conditional_t<K % C == 0, std::true_type, std::false_type> {};

template <typename T>
constexpr T default_cost();

template <>
constexpr size_t default_cost() { return std::numeric_limits<size_t>::max() / 2; }

template <typename CompType, typename CostType = size_t>
class sssp_solver {
 public:
  using cost_type = CostType;
  using comp_type = CompType;

  struct Vertex {
    cost_type cost;
    size_t len;
    size_t lz_ofs;
    comp_type type;
  };

  using vertex_type = Vertex;
  static constexpr cost_type default_cost = default_cost<cost_type>();

  sssp_solver() {}
  sssp_solver(const size_t n) : vertex(n + 1) {
    for (size_t i = 1; i <= n; ++i) (*this)[i].cost = default_cost;
    (*this)[0].cost = cost_type(0);
  }

  template <typename Func>
  void update_lz_table(size_t adr, std::span<const size_t> table, encode::lz_data lz, Func func, comp_type ty) {
    const cost_type base_cost = (*this)[adr].cost;
    for (size_t i = 0; i < table.size(); ++i) {
      const size_t l = table[i];
      if (l > lz.len || adr + l >= size()) break;
      const cost_type curr_cost = base_cost + func(i);
      if (curr_cost >= (*this)[adr + l].cost) continue;
      (*this)[adr + l].cost = curr_cost;
      (*this)[adr + l].len = l;
      (*this)[adr + l].lz_ofs = lz.ofs;
      (*this)[adr + l].type = ty;
    }
  }

  template <typename Pred = std::greater_equal<cost_type>, typename Func>
  void update(size_t adr, size_t fr, size_t to, Func func,
      comp_type ty, cost_type base_cost = default_cost, size_t optional = 0) {
    constexpr auto skip = Pred();
    to = std::min(to, size() > adr ? (size() - 1) - adr : 0);
    assert(fr > 0);
    if (base_cost == default_cost) base_cost = (*this)[adr].cost;
    for (size_t i = to; i >= fr; --i) {
      cost_type curr_cost = base_cost + func(i);
      if (skip(curr_cost, (*this)[adr + i].cost)) {
        if ((is_constant<Func>::value || is_linear<Func>::value) &&
             (*this)[adr + i].type == ty) {
          break;
        }
        continue;
      }
      (*this)[adr + i].cost = curr_cost;
      (*this)[adr + i].lz_ofs = optional;
      (*this)[adr + i].len = i;
      (*this)[adr + i].type = ty;
    }
  }

  template <typename Pred = std::greater_equal<cost_type>, typename Func>
  void update(size_t adr, size_t fr, size_t to, size_t len, Func func,
      comp_type ty, cost_type base_cost = default_cost, size_t optional = 0) {
    update<Pred, Func>(adr, fr, std::min(to, len), func, ty, base_cost, optional);
  }

  template <typename Pred = std::greater_equal<cost_type>, typename Func>
  void update_lz(size_t adr, size_t fr, size_t to, encode::lz_data lz, Func func,
      comp_type ty, cost_type base_cost = default_cost) {
    to = std::min(to, lz.len);
    update<Pred, Func>(adr, fr, std::min(to, lz.len), func, ty, base_cost, lz.ofs);
  }

  template <size_t K, typename Func>
  void update_k(size_t adr, size_t fr, size_t to, Func func, comp_type ty, size_t optional = 0) {
    to = std::min(to, size() > adr ? (size() - 1) - adr : 0);
    if (to < fr) return;
    to = fr + (to - fr) / K * K;
    cost_type base_cost = (*this)[adr].cost;
    for (size_t i = to; ptrdiff_t(i) >= ptrdiff_t(fr); i -= K) {
      cost_type curr_cost = base_cost + func(i);
      if (curr_cost >= (*this)[adr + i].cost) {
        if ((is_constant<Func>::value || is_linear_k<Func, K>::value) &&
             (*this)[adr + i].type == ty) {
          break;
        }
        continue;
      }
      (*this)[adr + i].cost = curr_cost;
      (*this)[adr + i].len = i;
      (*this)[adr + i].lz_ofs = optional;
      (*this)[adr + i].type = ty;
    }
  }

  template <size_t K, typename Func>
  void update_k(size_t adr, size_t fr, size_t to,
      size_t max_len, Func func, comp_type ty, size_t optional = 0) {
    update_k<K>(adr, fr, std::min(max_len, to), func, ty, optional);
  }

  cost_type total_cost() const {
    return vertex.back().cost;
  }

  std::vector<vertex_type> commands(ptrdiff_t start=0) const {
    std::vector<vertex_type> ret;
    ptrdiff_t adr = size() - 1;
    while (adr > start) {
      auto cmd = (*this)[adr];
      assert(cmd.len > 0);
      adr -= cmd.len;
      ret.emplace_back(cmd);
    }
    assert(adr == start);
    std::reverse(ret.begin(), ret.end());
    return ret;
  }

  const vertex_type& operator [] (size_t i) const {
    return vertex[i];
  }

  vertex_type& operator [] (size_t i) {
    return vertex[i];
  }

  size_t size() const {
    return vertex.size();
  }

 private:
  std::vector<vertex_type> vertex;
};

} // namespace sfc_comp
