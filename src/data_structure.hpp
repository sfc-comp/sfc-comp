#pragma once

#include <algorithm>
#include <vector>

#include <bit>

namespace sfc_comp {

template <typename ElementType, typename IndexType = uint32_t>
requires std::integral<ElementType> && std::unsigned_integral<IndexType> &&
         (sizeof(ElementType) <= sizeof(IndexType))
class suffix_array {
public:
  using element_type = ElementType;
  using index_type = IndexType;
  using signed_index_type = std::make_signed_t<index_type>;

  suffix_array() = default;
  suffix_array(std::span<const element_type> input) : input(input), sa(input.size()) {
    const size_t n = input.size();
    if (n == 0) return;

    for (size_t i = 0; i < n; ++i) sa[i] = i;
    std::vector<signed_index_type> inv_sa(n), ninv_sa(n);
    std::vector<index_type> psa(n);
    for (size_t i = 0; i < n; ++i) inv_sa[i] = input[i];

    std::sort(sa.begin(), sa.end(), [&] (const index_type a, const index_type b) {
      return inv_sa[a] < inv_sa[b] || (inv_sa[a] == inv_sa[b] && a > b);
    });

    std::vector<signed_index_type>& offsets = inv_sa;
    for (size_t l = 1; l < n; l <<= 1) {
      size_t lh = l >> 1;
      index_type ps = sa[0]; signed_index_type u = inv_sa[ps];
      ninv_sa[ps] = 0;
      size_t r = 1;
      for (size_t i = 1; i < n; ++i) {
        index_type cs = sa[i]; signed_index_type v = inv_sa[cs];
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
        signed_index_type k = psa[i] - l;
        if (k >= 0) sa[offsets[ninv_sa[k]]++] = k;
      }
      std::swap(inv_sa, ninv_sa);
    }
  }

  std::pair<std::vector<index_type>, std::vector<index_type>> lcp_rank() const {
    const size_t n = sa.size();
    std::vector<index_type> lcp(n, 0);
    std::vector<index_type> isa(n);
    for (size_t i = 0; i < n; ++i) isa[sa[i]] = i;
    size_t h = 0;
    for (size_t i = 0; i < n; ++i) {
      size_t r = isa[i];
      if (r + 1 >= n) continue;
      size_t j = sa[r + 1];
      for (size_t o = std::max(i, j); o + h < n && input[i + h] == input[j + h]; ++h);
      lcp[r] = h;
      if (h > 0) --h;
    }
    return std::make_pair(std::move(lcp), std::move(isa));
  }

  constexpr size_t size() const { return sa.size(); }

  index_type operator [] (size_t i) const {
    return sa[i];
  }

private:
  std::span<const element_type> input;
  std::vector<index_type> sa;
};

template <typename T>
concept monoid = requires (const typename T::value_type& a, const typename T::value_type& b) {
  typename T::value_type;
  { T::iden() } -> std::convertible_to<typename T::value_type>;
  { T::op(a, b) } -> std::convertible_to<typename T::value_type>;
};

template <typename T>
requires monoid<T>
class segment_tree {
public:
  using value_type = typename T::value_type;
  static constexpr value_type iden = T::iden();

  segment_tree() = default;
  segment_tree(const size_t n)
      : n(n), n2(n == 0 ? 0 : std::bit_ceil(n)), tree(2 * n2, T::iden()) {}

  segment_tree(std::span<const value_type> arr) : segment_tree(arr.size()) {
    init([&](size_t i) { return arr[i]; });
  }

  template <typename Func>
  requires std::convertible_to<std::invoke_result_t<Func, size_t>, value_type>
  void init(Func&& func) {
    if (n == 0) return;
    for (size_t i = 0; i < n; ++i) tree[n2 + i] = func(i);
    for (size_t i = n2 - 1; i > 0; --i) tree[i] = T::op(tree[2 * i], tree[2 * i + 1]);
  }

  template <typename Func>
  requires std::predicate<Func, value_type>
  std::pair<size_t, size_t> find_range(size_t i, Func&& func) const {
    size_t right = find_right(i, func);
    size_t left = find_left(i, func);
    return std::make_pair(left + 1, right);
  }

  template <typename Func>
  requires std::predicate<Func, value_type>
  ptrdiff_t find_right(size_t i, Func&& func) const {
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
  requires std::predicate<Func, value_type>
  ptrdiff_t find_left(size_t i, Func&& func) const {
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

  void reset(size_t k) {
    return update(k, T::iden());
  }

  value_type fold(size_t lo, size_t hi) const {
    value_type ret_l = T::iden(), ret_r = T::iden();
    for (lo += n2, hi += n2; lo < hi; lo >>= 1, hi >>= 1) {
      if (lo & 1) ret_l = T::op(ret_l, tree[lo++]);
      if (hi & 1) ret_r = T::op(tree[--hi], ret_r);
    }
    return T::op(ret_l, ret_r);
  }

  size_t size() const {
    return n;
  }

  std::span<const value_type> nodes() const {
    return std::span(tree.data(), tree.size());
  }

private:
  size_t n;
  size_t n2;
  std::vector<value_type> tree;
};

template <typename T>
requires std::unsigned_integral<T>
class wavelet_matrix {
 public:
  using value_type = T;
  static constexpr value_type nval = value_type(-1);

 private:
  using block_type = uint64_t;
  static constexpr size_t block_size = sizeof(block_type);

  struct counter {
    size_t cumu;
    block_type f;
  };

  size_t rank(size_t b, size_t i) const {
    const auto c = bit_vectors[b * vector_size + (i / block_size)];
    const auto mask = (block_type(1) << (i % block_size)) - 1;
    return c.cumu + std::popcount(c.f & mask);
  }

 public:
  wavelet_matrix() = default;
  wavelet_matrix(std::span<const value_type> input)
      : n(input.size()), vector_size(1 + n / block_size) {
    max_v = 0;
    for (const auto v : input) max_v = std::max(max_v, v);
    bit_width = std::bit_width(max_v);
    bit_vectors.resize(vector_size * bit_width);
    zeros.resize(bit_width);
    build(input);
  }

  value_type kth(size_t beg, size_t end, size_t k) const {
    if (end - beg <= k) return nval;
    value_type ret = 0;
    for (size_t b = bit_width; b-- > 0; ) {
      const size_t rb = rank(b, beg), re = rank(b, end), c = re - rb;
      if (c <= k) {
        beg += zeros[b] - rb; end += zeros[b] - re;
        k -= c; ret |= value_type(1) << b;
      } else {
        beg = rb; end = re;
      }
    }
    return ret;
  }

  size_t count_lt(size_t beg, size_t end, value_type v) const {
    if (v > max_v) return end - beg;
    if (v == 0) return 0;
    size_t ret = 0;
    for (size_t b = bit_width; beg < end && b-- > 0; ) {
      const size_t rb = rank(b, beg), re = rank(b, end);
      const auto mask = value_type(1) << b;
      if (v & mask) {
        ret += re - rb; v -= mask;
        if (v == 0) break;
        beg += zeros[b] - rb; end += zeros[b] - re;
      } else {
        beg = rb; end = re;
      }
    }
    return ret;
  }

 private:
  void build(std::span<const value_type> input) {
    std::vector<value_type> s(input.begin(), input.end()), t(n);
    for (size_t b = bit_width; b-- > 0; ) {
      size_t si = 0, ti = 0;
      auto bv = std::span(&bit_vectors[b * vector_size], vector_size);
      for (size_t i = 0; i < input.size(); ++i) {
        const auto v = s[i];
        if ((v >> b) & 1) {
          t[ti++] = v;
        } else {
          s[si++] = v;
          bv[i / block_size].f |= block_type(1) << (i % block_size);
        }
      }
      bv[0].cumu = 0;
      for (size_t i = 1; i < vector_size; ++i) {
        bv[i].cumu = bv[i - 1].cumu + std::popcount(bv[i - 1].f);
      }
      zeros[b] = si;
      std::copy_n(t.begin(), ti, s.begin() + si);
    }
  }

 private:
  size_t n;
  size_t vector_size;
  value_type max_v;
  size_t bit_width;
  std::vector<counter> bit_vectors;
  std::vector<size_t> zeros;
};

template <typename T, size_t WindowSize, typename Compare = std::greater<T>>
class sliding_window_max {
 public:
  static constexpr size_t capacity = WindowSize;
  using value_type = T;

  struct pair {
    value_type val;
    size_t key;
  };

  sliding_window_max() : head(0), tail(0) {}
  void add(size_t i, T val) {
    while (head < tail && Compare()(val, que[tail - 1].val)) --tail;
    que[tail++] = {val, i};
  }

  pair get() const { return que[head]; }
  void pop(size_t i) { if (que[head].key == i) ++head; }
  bool empty() const { return head == tail; }

 private:
  size_t head;
  size_t tail;
  std::array<pair, WindowSize> que;
};

} // namespace sfc_comp
