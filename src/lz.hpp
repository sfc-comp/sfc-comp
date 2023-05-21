#pragma once

#include <limits>
#include <span>

#include "data_structure.hpp"

namespace sfc_comp {

template <typename T, T Iden = std::numeric_limits<T>::min()>
struct range_max {
  using value_type = T;
  static constexpr T iden() { return Iden; }
  static constexpr T op(const value_type& l, const value_type& r) { return std::max<T>(l, r); }
};

template <typename T, T Iden = std::numeric_limits<T>::max()>
struct range_min {
  using value_type = T;
  static constexpr T iden() { return Iden; }
  static constexpr T op(const value_type& l, const value_type& r) { return std::min<T>(l, r); }
};

namespace encode {

struct lz_data {
  constexpr bool operator == (const lz_data& rhs) const = default;
  constexpr bool operator > (const lz_data& rhs) const {
    return len > rhs.len || (len == rhs.len && ofs > rhs.ofs);
  }
  size_t ofs;
  size_t len;
};

namespace lz {

template <typename U, typename S = std::make_signed_t<U>>
requires std::integral<U>
encode::lz_data find_left(size_t adr, size_t i, size_t d, size_t min_len,
    std::span<const U> lcp_node, std::span<const S> ofs_node) {
  if (i == 0) return {};
  const size_t width = lcp_node.size() / 2;
  const auto found = [&](size_t k) { return ofs_node[k] + ptrdiff_t(d) >= ptrdiff_t(adr); };
  U lcp = std::numeric_limits<U>::max();
  const auto quit = [&](size_t k) { return lcp_node[k] < lcp && (lcp = lcp_node[k]) < min_len; };

  size_t lo = i - 1, hi = i, k = lo + width;
  while (lo > 0 && !found(k)) {
    if (quit(k)) return {};
    size_t diff = hi - lo;
    if (!(k & 1)) hi = lo, lo -= 2 * diff, k = (k >> 1) - 1;
    else lo -= diff, hi -= diff, --k;
  }
  if (lo == 0 && !found(k)) return {};
  while (k < width) {
    size_t mi = (lo + hi) >> 1;
    if (found(2 * k + 1)) lo = mi, k = 2 * k + 1;
    else {
      if (quit(2 * k + 1)) return {};
      hi = mi, k = k * 2;
    }
  }
  if (quit(k)) return {};
  return {size_t(ofs_node[lo + width]), lcp};
}

template <typename U, typename S = std::make_signed_t<U>>
requires std::integral<U>
encode::lz_data find_right(size_t adr, size_t i, size_t d, size_t min_len,
    std::span<const U> lcp_node, std::span<const S> ofs_node) {
  const size_t width = lcp_node.size() / 2;
  const auto found = [&](size_t k) { return ofs_node[k] + ptrdiff_t(d) >= ptrdiff_t(adr); };
  U lcp = std::numeric_limits<U>::max();
  const auto quit = [&](size_t k) { return lcp_node[k] < lcp && (lcp = lcp_node[k]) < min_len; };

  size_t lo = i, hi = i + 1, k = lo + width;
  while (hi < width && !found(k)) {
    if (quit(k)) return {};
    size_t diff = hi - lo;
    if (k & 1) lo = hi, hi += 2 * diff, k = (k + 1) >> 1;
    else hi += diff, lo += diff, ++k;
  }
  if (hi == width && !found(k)) return {};
  while (k < width) {
    size_t mi = (lo + hi) >> 1;
    if (found(2 * k)) hi = mi, k = 2 * k;
    else {
      if (quit(2 * k)) return {};
      lo = mi, k = 2 * k + 1;
    }
  }
  return {size_t(ofs_node[lo + width]), lcp};
}

template <typename U, typename S = std::make_signed_t<U>>
requires std::integral<U>
encode::lz_data find(size_t adr, size_t rank, size_t max_dist, size_t min_len,
    std::span<const U> lcp_node, std::span<const S> ofs_node) {
  const auto left = find_left(adr, rank, max_dist, min_len, lcp_node, ofs_node);
  const auto right = find_right(adr, rank, max_dist, min_len, lcp_node, ofs_node);
  return left.len >= right.len ? left : right; // [TODO] choose the close one.
}

template <typename U, typename S = std::make_signed_t<U>>
requires std::integral<U>
encode::lz_data find_closest(size_t adr, size_t rank, size_t max_dist, size_t min_len, size_t max_len,
    const segment_tree<range_min<U>>& lcp, const segment_tree<range_max<S>>& seg) {
  auto ret = find(adr, rank, max_dist, min_len, lcp.nodes(), seg.nodes());
  if (ret.len > 0) {
    if (ret.len > max_len) ret.len = max_len;
    const auto r = lcp.find_range(rank, [&](size_t len) { return len >= ret.len; });
    ret.ofs = seg.fold(r.first, r.second + 1);
  }
  return ret;
}

template <typename U, typename Elem>
requires std::integral<U>
encode::lz_data find(size_t i, size_t j, size_t rank, const wavelet_matrix<U>& wm,
    const segment_tree<range_min<U>>& lcp, const suffix_array<Elem, U>& sa) {
  const auto k = wm.count_lt(i, j, rank);
  encode::lz_data ret = {};
  if (k > 0) {
    const auto rank_l = wm.kth(i, j, k - 1);
    const auto len_l = lcp.fold(rank_l, rank);
    if (len_l > ret.len) ret = {sa[rank_l], len_l};
  }
  if (k < (j - i)) {
    const auto rank_r = wm.kth(i, j, k);
    const auto len_r = lcp.fold(rank, rank_r);
    // [TODO] choose the close one.
    if (len_r > ret.len) ret = {sa[rank_r], len_r};
  }
  return ret;
}

template <typename Func>
requires std::convertible_to<std::invoke_result_t<Func, size_t>, encode::lz_data>
encode::lz_data find_non_overlapping(size_t adr_l, size_t adr, Func&& find_lz, encode::lz_data prev) {
  constexpr auto overlapped = [](size_t i, const encode::lz_data& res) {
    return res.len > 0 && res.ofs + res.len > i;
  };
  if (prev.len >= 1) prev.len -= 1, prev.ofs += 1;
  auto ret = find_lz(adr - (std::max<size_t>(prev.len, 1) - 1));
  if (!overlapped(adr, ret)) return ret;
  size_t len_hi = std::min(adr - adr_l, ret.len);
  ret.len = adr - ret.ofs;
  while (ret.len < len_hi) {
    const size_t len = (ret.len + len_hi + 1) / 2;
    auto lz = find_lz(adr - (len - 1));
    if (overlapped(adr, lz)) lz.len = adr - lz.ofs;
    if (lz > ret) ret = lz;
    if (lz.len < len) len_hi = len - 1;
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

  lz_helper(std::span<const uint8_t> input, bool updated = false) : n(input.size()) {
    const auto sa = suffix_array<uint8_t>(input);
    const auto [lcp, rank]= sa.lcp_rank();
    this->rank = std::move(rank);
    this->lcp = decltype(this->lcp)(lcp);
    this->seg = decltype(seg)(n);
    if (updated) this->seg.init([&](size_t i) { return sa[i]; });
  }

  encode::lz_data find(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg.nodes());
  }

  encode::lz_data find_closest(size_t pos, size_t max_dist, size_t min_len, size_t max_len) const {
    return encode::lz::find_closest(pos, rank[pos], max_dist, min_len, max_len, lcp, seg);
  }

  void reset(size_t i) {
    seg.reset(rank[i]);
  }

  void add_element(size_t i) {
    seg.update(rank[i], i);
  }

private:
  const size_t n;
  std::vector<index_type> rank;
  segment_tree<range_max<signed_index_type>> seg;
  segment_tree<range_min<index_type>> lcp;
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
  lz_helper_c(std::span<const uint8_t> input, bool updated = false) : n(input.size()) {
    const auto in = complement_appended(input);
    const auto sa = suffix_array<int16_t>(in);
    const auto [lcp, rank] = sa.lcp_rank();
    this->rank = std::move(rank);
    this->lcp = decltype(this->lcp)(lcp);
    seg = decltype(seg)(rank.size());
    if (updated) {
      seg.init([&](size_t i) { return sa[i] < n ? sa[i] : seg.iden; });
    }
    seg_c = decltype(seg_c)(rank.size());
    if (updated) {
      seg_c.init([&](size_t i) { return sa[i] >= n + 1 ? sa[i] - (n + 1) : seg_c.iden; });
    }
  }

public:
  encode::lz_data find(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg.nodes());
  }

  encode::lz_data find_c(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg_c.nodes());
  }

  void add_element(size_t i) {
    seg.update(rank[i], signed_index_type(i));
    seg_c.update(rank[i + n + 1], signed_index_type(i));
  }

  void reset(size_t i) {
    seg.reset(rank[i]);
    seg_c.reset(rank[i + n + 1]);
  }

private:
  const size_t n;
  std::vector<index_type> rank;
  segment_tree<range_max<signed_index_type>> seg, seg_c;
  segment_tree<range_min<index_type>> lcp;
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
  std::vector<int16_t> hvflip_appended(std::span<const uint8_t> input) const {
    std::vector<int16_t> ret(3 * n + 2);
    for (size_t i = 0; i < n; ++i) ret[i] = input[i];
    ret[n] = -1;
    for (size_t i = 0; i < n; ++i) ret[i + n + 1] = bit_reversed[input[i]];
    ret[2 * n + 1] = -2;
    for (size_t i = 0; i < n; ++i) ret[i + 2 * n + 2] = input[n - 1 - i];
    return ret;
  }

public:
  lz_helper_kirby(std::span<const uint8_t> input, bool updated = false) : n(input.size()) {
    const auto in = hvflip_appended(input);
    const auto sa = suffix_array<int16_t>(in);
    const auto [lcp, rank] = sa.lcp_rank();
    this->rank = std::move(rank);
    this->lcp = decltype(this->lcp)(lcp);
    seg = decltype(seg)(rank.size());
    seg_h = decltype(seg_h)(rank.size());
    seg_v = decltype(seg_v)(rank.size());
    if (updated) {
      seg.init([&](size_t i) { return sa[i] < n ? sa[i] : seg.iden; });
      seg_h.init([&](size_t i) { return n + 1 <= sa[i] && sa[i] < 2 * n + 1 ? sa[i] - (n + 1) : seg_h.iden; });
      seg_v.init([&](size_t i) { return sa[i] >= 2 * n + 2 ? (3 * n + 1) - sa[i] : seg_v.iden; });
    }
  }

public:
  encode::lz_data find(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg.nodes());
  }

  encode::lz_data find_h(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg_h.nodes());
  }

  encode::lz_data find_v(size_t pos, size_t max_dist, size_t min_len) const {
    return encode::lz::find(pos, rank[pos], max_dist, min_len, lcp.nodes(), seg_v.nodes());
  }

  void add_element(size_t i) {
    seg.update(rank[i], i);
    seg_h.update(rank[i + n + 1], i);
    seg_v.update(rank[3 * n + 1 - i], i);
  }

  void reset(size_t i) {
    seg.reset(rank[i]);
    seg_h.reset(rank[i + n + 1]);
    seg_v.reset(rank[3 * n + 1 - i]);
  }

private:
  const size_t n;
  std::vector<index_type> rank;
  segment_tree<range_min<index_type>> lcp;
  segment_tree<range_max<signed_index_type>> seg, seg_h, seg_v;
};

template <typename U = uint32_t>
requires std::unsigned_integral<U>
class non_overlapping_lz_helper {
 public:
  using index_type = U;
  using signed_index_type = std::make_signed_t<index_type>;

  non_overlapping_lz_helper(std::span<const uint8_t> input) : n(input.size()), sa(input) {
    const auto [lcp, rank]= sa.lcp_rank();
    this->rank = std::move(rank);
    this->wm = decltype(wm)(this->rank);
    this->lcp = decltype(this->lcp)(lcp);
  }

  encode::lz_data find_non_overlapping(const size_t adr, const size_t max_dist,
      const encode::lz_data prev = {}) const {
    const size_t adr_l = (adr < max_dist) ? 0 : adr - max_dist;
    const size_t rank = this->rank[adr];
    return encode::lz::find_non_overlapping(adr_l, adr, [&](size_t adr_r) {
      return encode::lz::find(adr_l, adr_r, rank, wm, lcp, sa);
    }, prev);
  }

  encode::lz_data find(const size_t adr, const size_t max_dist) const {
    const size_t adr_l = (adr < max_dist) ? 0 : adr - max_dist;
    return encode::lz::find(adr_l, adr, rank[adr], wm, lcp, sa);
  }

 private:
  const size_t n;
  suffix_array<uint8_t> sa;
  std::vector<index_type> rank;
  wavelet_matrix<index_type> wm;
  segment_tree<range_min<index_type>> lcp;
};

struct vrange {
  size_t min;
  size_t max;
  size_t bitlen;
  uint64_t val;
  uint64_t mask = -1;
};

struct vrange_min {
  size_t min;
  size_t bitlen;
  uint64_t val;
  uint64_t mask = -1;
};

template <size_t N>
constexpr std::array<vrange, N> to_vranges(vrange_min (&&a)[N], size_t max_len) {
  return create_array<vrange, N>([&](size_t i) {
    return vrange(a[i].min, (i + 1 == N) ? max_len : a[i + 1].min - 1, a[i].bitlen, a[i].val, a[i].mask);
  });
}

namespace encode::lz {

template <typename MaxOffset, typename Func>
requires std::convertible_to<std::invoke_result_t<MaxOffset, size_t>, size_t> &&
         std::convertible_to<std::invoke_result_t<Func, size_t>, encode::lz_data>
void find_all(size_t i, size_t o_size, const size_t lz_min_len,
    std::span<encode::lz_data> dest, MaxOffset&& max_ofs, Func&& find_lz) {
  for (ptrdiff_t oi = o_size - 1; oi >= 0; ) {
    auto res_lz = find_lz(max_ofs(oi));
    if (res_lz.len < lz_min_len) res_lz = {0, 0};
    do {
      dest[oi--] = res_lz;
    } while (oi >= 0 && (res_lz.len < lz_min_len || (i - res_lz.ofs) <= max_ofs(oi)));
  }
}

template <typename Func>
void find_all(size_t i, std::span<const size_t> max_offsets, const size_t lz_min_len,
    std::span<encode::lz_data> dest, Func&& find_lz) {
  return find_all(i, max_offsets.size(), lz_min_len, dest,
                  [&](size_t oi) { return max_offsets[oi]; }, std::forward<Func>(find_lz));
}

template <typename Func>
void find_all(size_t i, std::span<const vrange> offsets, const size_t lz_min_len,
    std::span<encode::lz_data> dest, Func&& find_lz) {
  return find_all(i, offsets.size(), lz_min_len, dest,
                  [&](size_t oi) { return offsets[oi].max; }, std::forward<Func>(find_lz));
}

} // namespace encode::lz

namespace detail {

template <typename LzFunc, typename Func>
requires std::convertible_to<std::invoke_result_t<LzFunc, size_t>, encode::lz_data> &&
         std::invocable<Func, size_t, size_t, size_t, encode::lz_data>
void for_each_best_lz(size_t adr, encode::lz_data res_lz, const size_t lz_min_len,
    std::span<const vrange> offsets, LzFunc&& find_lz, Func&& update) {
  using encode::lz_data;
  ptrdiff_t best_oi = -1; lz_data best_lz = {0, 0};
  size_t best = std::numeric_limits<size_t>::max();
  for (ptrdiff_t oi = offsets.size() - 1; res_lz.len >= lz_min_len; ) {
    const size_t d = adr - res_lz.ofs;
    while (oi >= 0 && d < offsets[oi].min) --oi;
    if (oi < 0) break;
    if (offsets[oi].bitlen <= best) {
      best_oi = oi; best = offsets[best_oi].bitlen; best_lz = res_lz;
    }
    lz_data nres_lz = (oi == 0) ? lz_data(0, 0) : find_lz(oi - 1);
    update(best_oi, nres_lz.len + 1, res_lz.len, best_lz);
    if (--oi < 0) break;
    res_lz = nres_lz;
  }
}

} // namespace detail

} // namespace sfc_comp
