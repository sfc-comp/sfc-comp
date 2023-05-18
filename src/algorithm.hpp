#pragma once

#include <cstddef>
#include <cassert>

#include <vector>
#include <limits>

#include <span>

#include "encode.hpp"
#include "data_structure.hpp"
#include "lz.hpp"

namespace sfc_comp {

template <typename CostType>
struct cost_traits;

template <>
struct cost_traits<size_t> {
  static constexpr size_t infinity() { return std::numeric_limits<size_t>::max() / 2; }
  static constexpr size_t unspecified() { return std::numeric_limits<size_t>::max() - 1; }
};

template <typename U, typename V>
concept add_able = requires (U u, V v) {
  {u + v} -> std::convertible_to<U>;
};

template <typename Tag>
requires std::equality_comparable<Tag> && std::convertible_to<Tag, uint32_t>
struct tag_ol {
  tag_ol() = default;
  tag_ol(Tag tag, uint64_t oi, uint64_t li) : tag(tag), oi(oi), li(li) {}
  Tag tag : 32;
  uint64_t oi : 16;
  uint64_t li : 16;
};

template <typename Tag>
requires std::equality_comparable<Tag> && std::convertible_to<Tag, uint32_t>
struct tag_l {
  tag_l() = default;
  tag_l(Tag tag, uint64_t li) : tag(tag), li(li) {}
  Tag tag : 32;
  uint64_t li : 16;
};

template <typename Tag>
requires std::equality_comparable<Tag> && std::convertible_to<Tag, uint32_t>
struct tag_o {
  tag_o() = default;
  tag_o(Tag tag, uint64_t oi) : tag(tag), oi(oi) {}
  Tag tag : 32;
  uint64_t oi : 16;
};

template <size_t A, size_t B, size_t D = 1>
struct linear {
  constexpr size_t operator()(size_t i) const { return (A * i + B) / D; }
};

template <size_t C>
struct constant : linear<0, C> {};

template <size_t Numer, size_t Denom = 1,
  typename Compare = std::greater<size_t>, typename CostType = size_t>
requires (Denom > 0)
class cost_window {
 public:
  using cost_type = CostType;
  static constexpr cost_type infinite_cost = cost_traits<cost_type>::infinity();
  static constexpr size_t nlen = std::numeric_limits<size_t>::max();

  struct len_cost {
    size_t len;
    cost_type cost;
  };

 private:
  struct value {
    constexpr auto operator < (const value& rhs) const {
      return cost < rhs.cost || (cost == rhs.cost && Compare()(index, rhs.index));
    }
    cost_type cost;
    size_t index;
  };
  static constexpr value iden = value(infinite_cost, nlen);

 public:
  cost_window() = default;
  cost_window(size_t size, size_t window_size, size_t beg = -2)
      : n(size), mask(std::bit_ceil((std::min(n + 1, window_size) + Denom - 1) / Denom) - 1) {
    for (size_t k = 0; k < Denom; ++k) segs[k] = segment_tree<range_min<value, iden>>(mask + 1);
    if (beg == size_t(-2)) beg = n;
    if (beg <= n) update(beg, 0);
  }

  void update(size_t i, cost_type cost) {
    segs[i % Denom].update((i / Denom) & mask, {cost + (i / Denom) * Numer, i});
  }

  constexpr cost_type operator [](size_t i) const {
    return segs[i % Denom][(i / Denom) & mask].cost + (i / Denom) * Numer;
  }

  len_cost find(size_t i, size_t fr, size_t to) const {
    if ((fr += i) > n) return {nlen, infinite_cost};
    to = std::min(n, i + to);
    const size_t d = (to - fr) / Denom + 1;
    const auto& seg = segs[fr % Denom];
    fr = (fr / Denom) & mask; to = fr + d;
    const auto res = (to <= mask + 1) ? seg.fold(fr, to)
                                      : std::min(seg.fold(fr, mask + 1), seg.fold(0, to & mask));
    if (res.cost >= infinite_cost) return {nlen, infinite_cost};
    return {res.index - i, res.cost - (res.index / Denom - (res.index - i) / Denom) * Numer};
  }

 private:
  size_t n;
  size_t mask;
  std::array<segment_tree<range_min<value, iden>>, Denom> segs;
};

template <typename TagType, typename CostType = size_t>
requires add_able<CostType, size_t>
class solver {
 public:
  using cost_type = CostType;
  using tag_type = TagType;
  static constexpr cost_type infinite_cost = cost_traits<cost_type>::infinity();

 public:
  struct node {
    constexpr node() = default;
    constexpr node(cost_type cost) : cost(cost) {}
    size_t lz_ofs() const { return arg; }
    cost_type cost;
    size_t len;
    size_t arg;
    tag_type type;
  };

 public:
  template <size_t Numer, size_t Denom, typename Less = std::greater<size_t>, typename C = cost_type>
  requires (Denom > 0) && std::convertible_to<C, cost_type>
  struct cmin : public cost_window<Numer, Denom, Less, C> {
   protected:
     using cost_window<Numer, Denom, Less, C>::update;
   public:
    cmin() : cost_window<Numer, Denom, Less, C>() {}
    cmin(std::span<const node> nd, size_t max_len, size_t dest = -2)
      : cost_window<Numer, Denom, Less, C>(nd.size() - 1, max_len, -1), nd(nd) {
      if (dest == size_t(-2)) dest = nd.size() - 1;
      if (dest < nd.size()) update(dest);
    }
    void update(size_t i) { update(i, nd[i].cost); }
   private:
    std::span<const node> nd;
  };

 public:
  solver() = default;
  solver(size_t n, size_t dest = -2) : n(n), nodes(n + 1, infinite_cost) {
    if (dest == size_t(-2)) dest = n;
    if (dest <= n) nodes[dest] = cost_type(0);
  }

  template <size_t Numer, size_t Denom = 1, typename Less = std::greater<size_t>, typename C = cost_type>
  cmin<Numer, Denom, Less, C> c(size_t max_len, size_t dest = -2) const {
    return cmin<Numer, Denom, Less, C>(this->nodes, max_len, dest);
  }

  template <typename Pred = std::less<cost_type>>
  void update_c(size_t adr, size_t l, cost_type cost, tag_type tag, size_t arg = 0) {
    if (auto& v = nodes[adr]; Pred()(cost, v.cost)) {
      v.cost = cost; v.len = l; v.arg = arg; v.type = tag;
    }
  }

  template <typename Pred = std::less<cost_type>>
  void update(size_t adr, size_t l, size_t c, tag_type tag, size_t arg = 0) {
    if (adr + l > n) return;
    update_c<Pred>(adr, l, nodes[adr + l].cost + c, tag, arg);
  }

  template <typename Pred = std::less<cost_type>, class RangeMin>
  void update(size_t adr, size_t fr, size_t to,
      const RangeMin& range_min, size_t c, tag_type tag, size_t arg = 0) {
    if (fr > to) return;
    auto res = range_min.find(adr, fr, to);
    if (res.len == range_min.nlen) return;
    update_c<Pred>(adr, res.len, res.cost + c, tag, arg);
  }

  template <typename Pred = std::less<cost_type>, class RangeMin>
  void update(size_t adr, size_t fr, size_t to, size_t len,
      const RangeMin& range_min, size_t c, tag_type tag, size_t arg = 0) {
    return update<Pred>(adr, fr, std::min(to, len), range_min, c, tag, arg);
  }

  template <typename Pred = std::less<cost_type>, class RangeMin>
  void update(size_t adr, size_t fr, size_t to, encode::lz_data lz,
      const RangeMin& range_min, size_t c, tag_type tag) {
    return update<Pred>(adr, fr, std::min(to, lz.len), range_min, c, tag, lz.ofs);
  }

  template <typename Pred = std::less<cost_type>, typename Cost>
  requires std::convertible_to<std::invoke_result_t<Cost, size_t>, size_t>
  void update_b(size_t adr, size_t fr, size_t to, Cost&& f, tag_type tag, size_t arg = 0) {
    for (size_t l = fr; l <= to; ++l) {
      if (adr + l > n) break;
      update_c<Pred>(adr, l, nodes[adr + l].cost + f(l), tag, arg);
    }
  }

  template <typename Pred = std::less<cost_type>, typename Cost>
  void update_b(size_t adr, size_t fr, size_t to, size_t len,
      Cost&& f, tag_type tag, size_t arg = 0) {
    return update_b<Pred>(adr, fr, std::min(len, to), std::forward<Cost>(f), tag, arg);
  }

  template <typename Pred = std::less<cost_type>, typename Cost>
  void update_b(size_t adr, size_t fr, size_t to, encode::lz_data lz,
      Cost&& f, tag_type tag) {
    return update_b<Pred>(adr, fr, std::min(lz.len, to), std::forward<Cost>(f), tag, lz.ofs);
  }

  template <typename Pred = std::less<cost_type>, class RangeMin, typename TagFunc>
  requires std::convertible_to<std::invoke_result_t<TagFunc, size_t>, tag_type>
  void update(size_t adr, std::span<const vrange> lranges, size_t len,
      const RangeMin& range_min, size_t c, TagFunc&& tag, size_t arg = 0) {
    for (size_t li = 0; li < lranges.size(); ++li) {
      const auto& l = lranges[li];
      if (len < l.min && adr + l.min > n) break;
      update<Pred>(adr, l.min, std::min(l.max, len), range_min, l.bitlen + c, tag(li), arg);
    }
  }

  template <typename Pred = std::less<cost_type>, class RangeMin, typename TagFunc>
  requires std::convertible_to<std::invoke_result_t<TagFunc, size_t>, tag_type>
  void update(size_t adr, std::span<const vrange> lranges,
      const RangeMin& range_min, size_t c, TagFunc&& tag, size_t arg = 0) {
    for (size_t li = 0; li < lranges.size(); ++li) {
      const auto& l = lranges[li];
      if (adr + l.min > n) break;
      update<Pred>(adr, l.min, l.max, range_min, l.bitlen + c, tag(li), arg);
    }
  }

  template <typename Pred = std::less<cost_type>, class RangeMin, typename TagFunc>
  requires std::convertible_to<std::invoke_result_t<TagFunc, size_t>, tag_type>
  void update(size_t adr, std::span<const vrange> lranges, encode::lz_data lz,
      const RangeMin& range_min, size_t c, TagFunc&& tag) {
    return update<Pred>(adr, lranges, lz.len, range_min, c, std::forward<TagFunc>(tag), lz.ofs);
  }

  template <typename Pred = std::less<cost_type>, typename Cost>
  requires std::convertible_to<std::invoke_result_t<Cost, size_t>, size_t>
  void update(size_t adr, std::span<const size_t> lens, encode::lz_data lz, Cost&& f, tag_type tag) {
    for (size_t li = 0; li < lens.size(); ++li) {
      const auto l = lens[li];
      if (adr + l > n || l > lz.len) break;
      update<Pred>(adr, l, f(li), tag, lz.ofs);
    }
  }

  template <typename Pred = std::less<cost_type>, class RangeMin, typename LzFunc, typename TagFunc>
  requires std::convertible_to<std::invoke_result_t<TagFunc, size_t, size_t>, tag_type>
  void update_matrix(size_t adr, std::span<const vrange> offsets, std::span<const vrange> lens,
      const RangeMin& range_min, size_t c, LzFunc&& find_lz, TagFunc&& tag) {
    if (lens.empty() || offsets.empty()) return;
    ptrdiff_t li = lens.size() - 1;
    const auto f = [&](size_t oi, size_t min_len, size_t max_len, encode::lz_data best_lz) {
      for (; li >= 0 && max_len < lens[li].min; --li);
      for (; li >= 0 && min_len <= lens[li].max; --li) {
        const auto& l = lens[li];
        update<Pred>(adr, std::max(min_len, l.min), std::min(max_len, l.max), best_lz, range_min,
               l.bitlen + offsets[oi].bitlen + c, tag(oi, li));
        if (min_len > l.min) break;
      }
    };
    return detail::for_each_best_lz(adr, find_lz(offsets.size() - 1), lens[0].min, offsets,
                                    std::forward<LzFunc>(find_lz), std::forward<decltype(f)>(f));
  }

  template <typename Pred = std::less<cost_type>, typename LzFunc, typename TagFunc, typename LenCostFunc>
  requires add_able<cost_type, std::invoke_result_t<LenCostFunc, size_t>> &&
           std::convertible_to<std::invoke_result_t<TagFunc, size_t, size_t>, tag_type>
  void update_matrix(size_t adr, std::span<const vrange> offsets, std::span<const size_t> lens,
      LenCostFunc&& len_cost, LzFunc&& find_lz, TagFunc&& tag) {
    if (lens.empty() || offsets.empty()) return;
    encode::lz_data res_lz = find_lz(offsets.size() - 1);
    res_lz.len = std::min(res_lz.len, n - adr); // needed when compressing chunks.
    ptrdiff_t li = (std::ranges::upper_bound(lens, res_lz.len) - lens.begin()) - 1;
    const auto f = [&](size_t oi, size_t min_len, size_t max_len, encode::lz_data best_lz) {
      for (; li >= 0 && max_len < lens[li]; --li);
      for (; li >= 0 && min_len <= lens[li]; --li) {
        update<Pred>(adr, lens[li], offsets[oi].bitlen + len_cost(li), tag(oi, li), best_lz.ofs);
      }
    };
    return detail::for_each_best_lz(adr, res_lz, lens[0], offsets,
                                    std::forward<LzFunc>(find_lz), std::forward<decltype(f)>(f));
  }

  const node& operator [] (size_t i) const { return nodes[i]; }

  struct path {
    template <typename Node>
    struct iterator {
      using iterator_category = std::forward_iterator_tag;
      using reference = Node&;
      using value_type = Node;
      using pointer = Node*;

      iterator(pointer p) : p(p) {}
      reference operator* () { return *p; }
      iterator& operator++ () {
        if (p->len == 0) throw std::logic_error("cmd.len == 0.");
        p += p->len;
        return *this;
      }
      iterator operator++ (int) {
        iterator ret = *this;
        ++(*this);
        return ret;
      }
      friend bool operator == (const iterator& lhs, const iterator& rhs) {
        return lhs.p == rhs.p;
      }
    private:
      pointer p;
    };
    using const_iterator = iterator<const node>;

    path(std::span<const node> nodes, size_t begin, size_t end)
        : b(&nodes[begin]), e(&nodes[end]) {
      if (begin > nodes.size() || end > nodes.size()) throw std::logic_error("invalid path.");
    }

    size_t size() const {
      size_t ret = 0;
      for (auto it = begin(); it != end(); ++it) ret += 1;
      return ret;
    }

    const_iterator begin() const { return b; }
    const_iterator end() const { return e; }

   private:
    const_iterator b, e;
  };

  path optimal_path(size_t begin = 0) const {
    return path(nodes, begin, n);
  }

  cost_type optimal_cost(size_t adr = 0) const {
    return nodes[adr].cost;
  }

 private:
  size_t n;
  std::vector<node> nodes;
};

} // namespace sfc_comp
