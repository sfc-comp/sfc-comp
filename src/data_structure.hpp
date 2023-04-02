#pragma once

#include <algorithm>
#include <vector>

#include "utility.hpp"

namespace sfc_comp {

template <typename T>
class suffix_array {
public:
  suffix_array() {}
  suffix_array(std::span<const T> input)
      : sa(input.size()) {

    const size_t n = input.size();
    if (n == 0) return;

    for (size_t i = 0; i < n; ++i) sa[i] = i;
    std::vector<int> inv_sa(n), ninv_sa(n), psa(n);
    for (size_t i = 0; i < n; ++i) inv_sa[i] = input[i];

    std::sort(sa.begin(), sa.end(), [&] (const int a, const int b) {
      return inv_sa[a] < inv_sa[b] || (inv_sa[a] == inv_sa[b] && a > b);
    });

    std::vector<int>& offsets = inv_sa;
    for (size_t l = 1; l < n; l <<= 1) {
      size_t lh = l >> 1;
      int ps = sa[0], u = inv_sa[ps];
      ninv_sa[ps] = 0;
      size_t r = 1;
      for (size_t i = 1; i < n; ++i) {
        int cs = sa[i], v = inv_sa[cs];
        if (u == v && ps + l < n && inv_sa[ps + lh] == inv_sa[cs + lh]) {
          ninv_sa[cs] = ninv_sa[ps];
        } else {
          ninv_sa[cs] = i; ++r;
        }
        ps = cs; u = v;
      }
      if (r == n) break;
      for (size_t i = 0; i < n; ++i) psa[i] = sa[i];
      for (size_t i = 0; i < n; ++i) offsets[i] = i;
      for (size_t i = 0; i < n; ++i) {
        int k = psa[i] - l;
        if (k >= 0) sa[offsets[ninv_sa[k]]++] = k;
      }
      std::swap(inv_sa, ninv_sa);
    }
  }

  std::pair<std::vector<int>, std::vector<int>> compute_lcp_and_rank(std::span<const T> input) const {
    size_t n = sa.size();
    std::vector<int> lcp(n, 0);
    std::vector<int> isa(n);
    for (size_t i = 0; i < n; ++i) {
      isa[sa[i]] = i;
    }
    size_t h = 0;
    for (size_t i = 0; i < n; ++i) {
      size_t r = isa[i];
      if (r + 1 >= n) continue;
      size_t j = sa[r + 1];
      for (size_t o = std::max(i, j); o + h < n && input[i + h] == input[j + h]; ++h);
      lcp[r] = h;
      if (h > 0) --h;
    }
    return std::make_pair(lcp, isa);
  }

  int operator [] (size_t i) const {
    return sa[i];
  }

private:
  std::vector<int> sa;
};

template <typename T>
class segment_tree {
public:
  using value_type = typename T::value_type;

  segment_tree() {}

  segment_tree(const size_t n)
      : n(n), n2(1 << ilog2(2 * n + 1)), tree(2 * n2, T::unit()) {}

  segment_tree(std::span<const value_type> arr)
      : n(arr.size()), n2(1 << ilog2(2 * n + 1)), tree(2 * n2, T::unit()) {
    for (size_t i = 0; i < n; ++i) tree[n2 + i] = arr[i];
    for (size_t i = n2 - 1; i > 0; --i) tree[i] = T::op(tree[2 * i], tree[2 * i + 1]);
  }

  template <typename Func>
  std::pair<size_t, size_t> find_range(size_t i, Func func) const {
    auto right = find_right(i, func);
    auto left = find_left(i, func);
    return std::make_pair(left + 1, right);
  }

  template <typename Func>
  ptrdiff_t find_right(size_t i, Func func) const {
    size_t lo = i, hi = i + 1, k = lo + n2;
    while (hi < n2 && func(tree[k])) {
      size_t diff = hi - lo;
      if (k & 1) lo = hi, hi += 2 * diff, k = (k + 1) >> 1;
      else hi += diff, lo += diff, ++k;
    }
    if (hi == n2 && func(tree[k])) return n;
    while (k < n2) {
      size_t mi = (lo + hi) >> 1;
      if (!func(tree[2 * k])) hi = mi, k = 2 * k;
      else lo = mi, k = 2 * k + 1;
    }
    return lo;
  }

  template <typename Func>
  ptrdiff_t find_left(size_t i, Func func) const {
    if (i == 0) return -1;
    size_t lo = i - 1, hi = i, k = lo + n2;
    while (lo > 0 && func(tree[k])) {
      size_t diff = hi - lo;
      if (!(k & 1)) hi = lo, lo -= 2 * diff, k = (k >> 1) - 1;
      else lo -= diff, hi -= diff, --k;
    }
    if (lo == 0 && func(tree[k])) return -1;
    while (k < n2) {
      size_t mi = (lo + hi) >> 1;
      if (!func(tree[2 * k + 1])) lo = mi, k = 2 * k + 1;
      else hi = mi, k = 2 * k;
    }
    return lo;
  }

  value_type operator [] (const size_t i) const {
    return tree[n2 + i];
  }

  void update(size_t k, value_type v) {
    k += n2;
    tree[k] = v;
    for (k >>= 1; k > 0; k >>= 1) tree[k] = T::op(tree[2 * k], tree[2 * k + 1]);
  }

  value_type fold(size_t lo, size_t hi) const {
    value_type ret_l = T::unit(), ret_r = T::unit();
    for (lo += n2, hi += n2; lo < hi; lo >>= 1, hi >>= 1) {
      if (lo & 1) ret_l = T::op(ret_l, tree[lo++]);
      if (hi & 1) ret_r = T::op(tree[--hi], ret_r);
    }
    return T::op(ret_l, ret_r);
  }

  size_t size() const {
    return n;
  }

private:
  size_t n;
  size_t n2;
  std::vector<value_type> tree;
};

} // namespace sfc_comp
