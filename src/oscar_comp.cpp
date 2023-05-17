#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> oscar_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = create_array<vrange, 15>([](size_t k) {
    return vrange(1 << k, (2 << k) - 1, 4 + k, (k + 1) << k);
  });
  static constexpr auto ulen_tab = create_array<vrange, 8>([](size_t k) {
    return (k == 0) ? vrange(1, 1, 3, 0)
                    : vrange((1 << (k - 1)) + 1, (1 << k), 3 + (k - 1), k << (k - 1));
  });
  static constexpr size_t min_len = 3;
  static constexpr auto len_tab = create_array<vrange, ulen_tab.size()>([](size_t k) {
    return (k == 0) ? vrange(min_len, min_len, 3, 0)
                    : vrange((1 << (k - 1)) + min_len, (1 << k) + (min_len - 1), 3 + (k - 1), k << (k - 1));
  });

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(len_tab.back().max);
  auto c8 = dp.c<8>(ulen_tab.back().max);

  for (size_t i = input.size(); i-- > 0; ) {
    lz_helper.reset(i);
    dp.update(i, ulen_tab, c8, 4, [&](size_t li) -> tag { return {uncomp, 0, li}; });
    dp.update_matrix(i, ofs_tab, len_tab, c0, 0,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i); c8.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(4);
  writer raw;
  size_t adr = 0, counter = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    const size_t d = adr - cmd.lz_ofs();
    switch (tag) {
    case uncomp: {
      const auto& l = ulen_tab[li];
      ret.write<bnh, bnh>({4, 0}, {l.bitlen, l.val + (cmd.len - l.min)});
      raw.write<d8n>({cmd.len, &input[adr]});
    } break;
    case lz: {
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[oi];
      ret.write<bnh>({o.bitlen, o.val + (d - o.min)});
      ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
    ++counter;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 32 == raw.size() * 8 + ret.bit_length());
  write16b(ret.out, 0, ret.size() - 4);
  write16b(ret.out, 2, counter);
  ret.extend(raw);
  return ret.out;
}

} // namespace sfc_comp
