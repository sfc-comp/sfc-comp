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

template <typename T, T Iden = std::numeric_limits<T>::min()>
struct range_max {
  using value_type = T;
  static T iden() { return Iden; }
  static T op(const value_type& l, const value_type& r) { return std::max<T>(l, r); }
};

template <typename T, T Iden = std::numeric_limits<T>::max()>
struct range_min {
  using value_type = T;
  static T iden() { return Iden; }
  static T op(const value_type& l, const value_type& r) { return std::min<T>(l, r); }
};

namespace encode {

struct lz_data {
  ptrdiff_t ofs;
  size_t len;
};

namespace lz {

template <typename U, typename S>
requires std::integral<U> && std::signed_integral<S>
encode::lz_data find_best(
    size_t adr, size_t index, size_t block_size,
    const segment_tree<range_min<U>>& seg_lcp,
    const segment_tree<range_max<S>>& seg) {
  auto r = seg.find_range(index, [&](S v) {
    return v + static_cast<S>(block_size) < static_cast<S>(adr);
  });
  U len = 0; S offset = -1;
  if (r.first > 0) {
    const U len_l = seg_lcp.fold(r.first - 1, index);
    if (len_l > len) len = len_l, offset = seg[r.first - 1];
  }
  if (r.second < seg.size()) {
    const U len_r = seg_lcp.fold(index, r.second);
    if (len_r > len) len = len_r, offset = seg[r.second];
  }
  return {offset, len};
}

template <typename U, typename S>
requires std::integral<U> && std::signed_integral<S>
encode::lz_data find_best_closest(
    size_t pos, size_t index, size_t block_size, size_t max_len,
    const segment_tree<range_min<U>>& seg_lcp,
    const segment_tree<range_max<S>>& seg) {
  auto ret = find_best(pos, index, block_size, seg_lcp, seg);
  if (ret.len > 0) {
    if (ret.len > max_len) ret.len = max_len;
    const auto r2 = seg_lcp.find_range(index, [&](size_t v) { return v >= ret.len; });
    ret.ofs = seg.fold(r2.first, r2.second + 1);
  }
  return ret;
}

} // namespace lz

} // namespace encode

template <typename U = uint32_t>
requires std::unsigned_integral<U>
class lz_helper {
public:
  using index_type = U;
  using signed_index_type = std::make_signed_t<index_type>;

  lz_helper(std::span<const uint8_t> input) : n(input.size()), segtree(n) {
    const auto [lcp, rank]= suffix_array<uint8_t>(input).lcp_rank();
    this->rank = std::move(rank);
    segtree_lcp = decltype(segtree_lcp)(lcp);
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
  std::vector<index_type> rank;
  segment_tree<range_max<signed_index_type>> segtree;
  segment_tree<range_min<index_type>> segtree_lcp;
};

template <typename U = uint32_t>
requires std::unsigned_integral<U>
class lz_helper_c {
public:
  using index_type = U;
  using signed_index_type = std::make_signed_t<index_type>;

private:
  std::vector<int16_t> complement_appended(std::span<const uint8_t> input) const {
    std::vector<int16_t> input_xor(2 * n + 1);
    for (size_t i = 0; i < n; ++i) input_xor[i] = input[i];
    input_xor[n] = -1;
    for (size_t i = 0; i < n; ++i) input_xor[i + n + 1] = input[i] ^ 0xff;
    return input_xor;
  }

public:
  lz_helper_c(std::span<const uint8_t> input) : n(input.size()) {
    const auto [lcp, rank] = suffix_array<int16_t>(complement_appended(input)).lcp_rank();
    this->rank = std::move(rank);
    segtree = decltype(segtree)(rank.size());
    segtree_c = decltype(segtree_c)(rank.size());
    segtree_lcp = decltype(segtree_lcp)(lcp);
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
    segtree.update(rank[i], signed_index_type(i));
    segtree_c.update(rank[i + n + 1], signed_index_type(i));
  }

private:
  const size_t n;
  std::vector<index_type> rank;
  segment_tree<range_max<signed_index_type>> segtree, segtree_c;
  segment_tree<range_min<index_type>> segtree_lcp;
};

constexpr auto bit_reversed = [] {
  std::array<uint8_t, 256> rev;
  for (size_t i = 0; i < rev.size(); ++i) {
    size_t b = i;
    b = (b & 0xf0) >> 4 | (b & 0x0f) << 4;
    b = (b & 0xcc) >> 2 | (b & 0x33) << 2;
    b = (b & 0xaa) >> 1 | (b & 0x55) << 1;
    rev[i] = b;
  }
  return rev;
}();

template <typename U = uint32_t>
requires std::unsigned_integral<U>
class lz_helper_kirby {
public:
  using index_type = U;
  using signed_index_type = std::make_signed_t<index_type>;

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
  lz_helper_kirby(std::span<const uint8_t> input) : n(input.size()) {
    const auto [lcp_h, rank_h] = suffix_array<int16_t>(hflip_appended(input)).lcp_rank();
    this->rank_h = std::move(rank_h);
    seg = decltype(seg)(rank_h.size());
    seg_h = decltype(seg_h)(rank_h.size());
    seg_lcp_h = decltype(seg_lcp_h)(lcp_h);

    const auto [lcp_v, rank_v] = suffix_array<int16_t>(vflip_appended(input)).lcp_rank();
    this->rank_v = std::move(rank_v);
    seg_v = decltype(seg_v)(rank_v.size());
    seg_lcp_v = decltype(seg_lcp_v)(lcp_v);
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
  std::vector<index_type> rank_h, rank_v;
  segment_tree<range_min<index_type>> seg_lcp_h, seg_lcp_v;
  segment_tree<range_max<signed_index_type>> seg, seg_h, seg_v;
};

template <size_t A, size_t B, size_t C>
struct LinearQ {
  constexpr size_t operator() (size_t i) const {
    return (A * i + B) / C;
  }
};

template <size_t A, size_t B>
struct Linear : LinearQ<A, B, 1> {};

template <size_t C>
struct Constant : Linear<0, C> {};

template <typename>
struct is_linear : std::false_type {};

template <size_t A, size_t B, size_t C>
struct is_linear<LinearQ<A, B, C>>
  : std::conditional_t<C == 1, std::true_type, std::false_type> {};

template <size_t A, size_t B>
struct is_linear<Linear<A, B>> : is_linear<LinearQ<A, B, 1>> {};

template <size_t N>
struct is_linear<Constant<N>> : is_linear<Linear<0, N>> {};

template <typename, size_t>
struct is_linear_k : std::false_type {};

template <size_t A, size_t B, size_t C, size_t K>
struct is_linear_k<LinearQ<A, B, C>, K>
  : std::conditional_t<K % C == 0, std::true_type, std::false_type> {};

template <size_t A, size_t B, size_t K>
struct is_linear_k<Linear<A, B>, K> : is_linear_k<LinearQ<A, B, 1>, K> {};

template <size_t N, size_t K>
struct is_linear_k<Constant<N>, K> : is_linear_k<Linear<0, N>, K> {};

template <typename CostType>
struct cost_traits;

template <>
struct cost_traits<size_t> {
  static constexpr size_t infinity() { return std::numeric_limits<size_t>::max() / 2; }
  static constexpr size_t unspecified() { return std::numeric_limits<size_t>::max(); }
};

template <typename CostType = size_t>
class uncomp_helper {
 public:
  using cost_type = CostType;
  static constexpr cost_type infinite_cost = cost_traits<cost_type>::infinity();
  static constexpr size_t nlen = std::numeric_limits<size_t>::max();

  struct len_cost {
    size_t len;
    cost_type cost;
  };

 private:
  struct indexed_cost {
    constexpr bool operator < (const indexed_cost& rhs) const {
      return cost < rhs.cost || (cost == rhs.cost && index < rhs.index);
    }
    cost_type cost;
    size_t index;
  };
  static constexpr indexed_cost iden = indexed_cost(infinite_cost, nlen);

 public:
  uncomp_helper(size_t size, size_t slope)
      : n(size), slope_(slope), tree_(size) {}

  void update(size_t i, cost_type cost) {
    tree_.update(i, {cost  + (n - i) * slope_, i});
  }

  void reset(size_t i) {
    tree_.update(i, iden);
  }

  void reset(size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) reset(i);
  }

  len_cost find(size_t i, size_t fr, size_t to) const {
    if (i < fr) return {nlen, infinite_cost};
    to = std::min(i, to);
    const auto res = tree_.fold(i - to, i - fr + 1);
    if (res.cost >= infinite_cost) return {nlen, infinite_cost};
    return {i - res.index, res.cost - (n - i) * slope_};
  }

 private:
  size_t n;
  size_t slope_;
  segment_tree<range_min<indexed_cost, iden>> tree_;
};

struct vrange {
  size_t min;
  size_t max;
  size_t bitlen;
  size_t val;
};

template <typename Tag>
requires std::equality_comparable<Tag> && std::convertible_to<Tag, uint16_t>
struct tag_ol {
  tag_ol() = default;
  tag_ol(Tag tag, uint64_t oi, uint64_t li) : tag(tag), oi(oi), li(li) {}
  constexpr bool operator == (const tag_ol& rhs) const {
    return tag == rhs.tag && li == rhs.li;
  }
  Tag tag : 16;
  uint32_t oi : 8;
  uint32_t li : 8;
};

template <typename Tag>
requires std::equality_comparable<Tag> && std::convertible_to<Tag, uint16_t>
struct tag_l {
  tag_l() = default;
  tag_l(Tag tag, uint64_t li) : tag(tag), li(li) {}
  constexpr bool operator == (const tag_l& rhs) const {
    return tag == rhs.tag && li == rhs.li;
  }
  Tag tag : 16;
  uint32_t li : 8;
};

template <typename Tag>
requires std::equality_comparable<Tag> && std::convertible_to<Tag, uint16_t>
struct tag_o {
  tag_o() = default;
  tag_o(Tag tag, uint64_t oi) : tag(tag), oi(oi) {}
  constexpr bool operator == (const tag_o& rhs) const {
    return tag == rhs.tag;
  }
  Tag tag : 16;
  uint32_t oi : 8;
};

template <typename U, typename V>
concept add_able = requires (U u, V v) {
  {u + v} -> std::convertible_to<U>;
};

template <typename TagType, typename CostType = size_t>
requires std::equality_comparable<TagType>
class sssp_solver {
 public:
  using cost_type = CostType;
  using tag_type = TagType;
  static constexpr cost_type infinite_cost = cost_traits<cost_type>::infinity();

 private:
  static constexpr cost_type unspecified = cost_traits<cost_type>::unspecified();

 public:
  struct Vertex {
    cost_type cost;
    size_t len;
    size_t lz_ofs;
    tag_type type;

    size_t val() const { return lz_ofs; }
  };

  using vertex_type = Vertex;

 public:
  sssp_solver() = default;
  sssp_solver(const size_t n, size_t begin = 0) : vertex(n + 1) {
    reset(0, n + 1);
    if (begin <= n) (*this)[begin].cost = cost_type(0);
  }

  void reset(size_t i) {
    (*this)[i].cost = infinite_cost;
  }

  void reset(size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) reset(i);
  }

  template <typename LzFunc, typename TagFunc>
  requires std::convertible_to<std::invoke_result_t<LzFunc, size_t>, encode::lz_data> &&
           std::convertible_to<std::invoke_result_t<TagFunc, size_t, size_t>, tag_type>
  void update_lz_matrix(size_t adr, std::span<const vrange> offsets, std::span<const vrange> lens,
      LzFunc&& find_lz, TagFunc&& tag, size_t c, cost_type base_cost = unspecified) {
    if (lens.empty() || offsets.empty()) return;
    const size_t lz_min_len = lens[0].min;
    if (base_cost == unspecified) base_cost = (*this)[adr].cost;

    using encode::lz_data;
    ptrdiff_t oi = offsets.size() - 1;
    ptrdiff_t li = lens.size() - 1;

    ptrdiff_t best_oi = -1;
    size_t best_bitlen = std::numeric_limits<size_t>::max();
    encode::lz_data best_lz = {0, 0};
    for (lz_data res_lz = find_lz(oi); ; ) {
      const size_t d = adr - res_lz.ofs;
      if (res_lz.len < lz_min_len) break;
      while (oi >= 0 && d < offsets[oi].min) --oi;
      if (oi < 0) break;
      if (offsets[oi].bitlen <= best_bitlen) {
        best_bitlen = offsets[oi].bitlen; best_oi = oi; best_lz = res_lz;
      }
      lz_data nres_lz = (oi == 0) ? lz_data(0, 0) : find_lz(oi - 1);
      for (; li >= 0 && res_lz.len < lens[li].min; --li);
      const size_t min_len = nres_lz.len + 1;
      for (; li >= 0 && min_len <= lens[li].max; --li) {
        const auto& l = lens[li];
        const auto cost = base_cost + (best_bitlen + l.bitlen + c);
        update_lz(adr, std::max(min_len, l.min), l.max, best_lz, Constant<0>(), tag(best_oi, li), cost);
        if (min_len > l.min) break;
      }
      if (--oi < 0) break;
      res_lz = std::move(nres_lz);
    }
  }

  template <typename LzFunc, typename TagFunc, typename LenCostFunc>
  requires std::convertible_to<std::invoke_result_t<LzFunc, size_t>, encode::lz_data> &&
           add_able<cost_type, std::invoke_result_t<LenCostFunc, size_t>> &&
           std::convertible_to<std::invoke_result_t<TagFunc, size_t, size_t>, tag_type>
  void update_lz_matrix(size_t adr, std::span<const vrange> offsets, std::span<const size_t> lens,
      LzFunc&& find_lz, LenCostFunc&& func, TagFunc&& tag, cost_type base_cost = unspecified) {
    if (lens.empty() || offsets.empty()) return;
    const size_t lz_min_len = lens[0];
    if (base_cost == unspecified) base_cost = (*this)[adr].cost;

    using encode::lz_data;
    ptrdiff_t oi = offsets.size() - 1;
    lz_data res_lz = find_lz(oi);
    ptrdiff_t li = (std::upper_bound(lens.begin(), lens.end(), res_lz.len) - lens.begin()) - 1;

    ptrdiff_t best_oi = -1;
    size_t best_bitlen = std::numeric_limits<size_t>::max();
    encode::lz_data best_lz = {0, 0};
    for (; ;) {
      const size_t d = adr - res_lz.ofs;
      if (res_lz.len < lz_min_len) break;
      while (oi >= 0 && d < offsets[oi].min) --oi;
      if (oi < 0) break;
      if (offsets[oi].bitlen <= best_bitlen) {
        best_bitlen = offsets[oi].bitlen; best_oi = oi; best_lz = res_lz;
      }
      lz_data nres_lz = (oi == 0) ? lz_data(0, 0) : find_lz(oi - 1);
      for (; li >= 0 && res_lz.len < lens[li]; --li);
      const size_t min_len = nres_lz.len + 1;
      for (; li >= 0 && min_len <= lens[li]; --li) {
        const size_t l = lens[li];
        const auto cost = base_cost + (best_bitlen + func(li));
        if (adr + l >= size()) continue;
        auto& target = (*this)[adr + l];
        if (cost >= target.cost) continue;
        target = {cost, l, size_t(best_lz.ofs), tag(best_oi, li)};
      }
      if (--oi < 0) break;
      res_lz = std::move(nres_lz);
    }
  }

  template <typename Func>
  requires add_able<cost_type, std::invoke_result_t<Func, size_t>>
  void update_lz_table(size_t adr, std::span<const size_t> table, encode::lz_data lz, Func&& func, tag_type tag) {
    const cost_type base_cost = (*this)[adr].cost;
    for (size_t i = 0; i < table.size(); ++i) {
      const size_t l = table[i];
      if (l > lz.len || adr + l >= size()) break;
      auto& target = (*this)[adr + l];
      const cost_type curr_cost = base_cost + func(i);
      if (curr_cost >= target.cost) continue;
      target.cost = curr_cost;
      target.len = l;
      target.lz_ofs = lz.ofs;
      target.type = tag;
    }
  }

  template <typename Skip = std::greater_equal<cost_type>>
  requires std::convertible_to<std::invoke_result_t<Skip, cost_type, cost_type>, bool>
  void update(size_t adr, size_t len, tag_type tag, cost_type cost, size_t arg = 0) {
    if (adr < len) return;
    constexpr auto skip = Skip();
    auto& target = (*this)[adr];
    if (skip(cost, target.cost)) return;
    target.cost = cost;
    target.len = len;
    target.lz_ofs = arg;
    target.type = tag;
  }

  void update_u(size_t adr, size_t len, tag_type tag, cost_type cost, size_t arg = 0) {
    return update<std::greater<cost_type>>(adr, len, tag, cost, arg);
  }

  template <typename Skip = std::greater_equal<cost_type>, typename Func>
  requires std::convertible_to<std::invoke_result_t<Skip, cost_type, cost_type>, bool> &&
           add_able<cost_type, std::invoke_result_t<Func, size_t>>
  void update(size_t adr, size_t fr, size_t to, Func&& func,
      tag_type tag, cost_type base_cost = unspecified, size_t arg = 0) {
    constexpr auto skip = Skip();
    to = std::min(to, size() > adr ? (size() - 1) - adr : 0);
    if (base_cost == unspecified) base_cost = (*this)[adr].cost;
    for (size_t i = to; ptrdiff_t(i) >= ptrdiff_t(fr); --i) {
      const cost_type curr_cost = base_cost + func(i);
      auto& target = (*this)[adr + i];
      if (skip(curr_cost, target.cost)) {
        if constexpr (is_linear<Func>::value) {
          if (target.type == tag) break;
        }
        continue;
      }
      target.cost = curr_cost;
      target.lz_ofs = arg;
      target.len = i;
      target.type = tag;
    }
  }

  template <typename Skip = std::greater_equal<cost_type>, typename Func>
  void update(size_t adr, size_t fr, size_t to, size_t len, Func&& func,
      tag_type tag, cost_type base_cost = unspecified, size_t arg = 0) {
    return update<Skip, Func>(adr, fr, std::min(to, len), std::forward<Func>(func), tag, base_cost, arg);
  }

  template <typename Skip = std::greater_equal<cost_type>, typename Func>
  void update_lz(size_t adr, size_t fr, size_t to, encode::lz_data lz, Func&& func,
      tag_type tag, cost_type base_cost = unspecified) {
    return update<Skip, Func>(adr, fr, std::min(to, lz.len), std::forward<Func>(func), tag, base_cost, lz.ofs);
  }

  template <size_t K, typename Func>
  requires (K > 0) &&
           add_able<cost_type, std::invoke_result_t<Func, size_t>>
  void update_k(size_t adr, size_t fr, size_t to, Func&& func, tag_type tag, size_t arg = 0) {
    to = std::min(to, size() > adr ? (size() - 1) - adr : 0);
    if (to < fr) return;
    to = fr + (to - fr) / K * K;
    const cost_type base_cost = (*this)[adr].cost;
    for (size_t i = to; ptrdiff_t(i) >= ptrdiff_t(fr); i -= K) {
      const cost_type curr_cost = base_cost + func(i);
      auto& target = (*this)[adr + i];
      if (curr_cost >= target.cost) {
        if constexpr (is_linear_k<Func, K>::value) {
          if (target.type == tag) break;
        }
        continue;
      }
      target.cost = curr_cost;
      target.len = i;
      target.lz_ofs = arg;
      target.type = tag;
    }
  }

  template <size_t K, typename Func>
  void update_k(size_t adr, size_t fr, size_t to,
      size_t max_len, Func&& func, tag_type tag, size_t arg = 0) {
    return update_k<K>(adr, fr, std::min(max_len, to), std::forward<Func>(func), tag, arg);
  }

  cost_type total_cost() const {
    return vertex.back().cost;
  }

  std::vector<vertex_type> commands(ptrdiff_t start=0) const {
    std::vector<vertex_type> ret;
    ptrdiff_t adr = size() - 1;
    while (adr > start) {
      auto cmd = (*this)[adr];
      if (cmd.len == 0) throw std::logic_error("cmd.len == 0");
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
