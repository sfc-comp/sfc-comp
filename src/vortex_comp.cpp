#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> vortex_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0x800000);

  enum method { none, uncomp, lz2, lz3, lz };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = to_vranges({
             // 00 (len == 2)
             // 01 (len == 3)
    {  4,  2, 0b10},
    {  5,  4, 0b110'0},
    {  7,  6, 0b1110'00},
    { 11, 12, 0b1111'00000000 + 0x0b}
  }, 255);

  static constexpr auto ofs_tab = to_vranges({
    {0x0001, 10, 0b11'00000000 + 0x0001},
    {0x0100, 14, 0b10'000000000000 + 0x0100},
    {0x1000, 17, 0b0'0000000000000000 + 0x1000}
  }, 0xffff);

  std::vector<uint8_t> input(in.rbegin(), in.rend());

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<tag> dp0(input.size()), dp1(input.size(), -1);

  for (size_t i = 0; i <= input.size(); ++i) {
    const auto cost0 = dp0[i].cost;
    if (cost0 < dp0.infinite_cost) {
      dp1.update(i, 0, 0, Constant<3>(), {none, 0, 0}, cost0);
    }
    if (i == input.size()) break;

    if (cost0 < dp0.infinite_cost) u_helper.update(i, cost0);
    auto u0 = u_helper.find(i + 1, 1, 6);
    dp1.update_u(i + 1, u0.len, {uncomp, 0, 0}, u0.cost + 3);
    auto u1 = u_helper.find(i + 1, 7, 0x16);
    dp1.update_u(i + 1, u1.len, {uncomp, 0, 1}, u1.cost + 8);
    auto u2 = u_helper.find(i + 1, 0x17, 0x0fff);
    dp1.update_u(i + 1, u2.len, {uncomp, 0, 2}, u2.cost + 14);

    const auto cost1 = dp1[i].cost;
    if (cost1 < dp1.infinite_cost) {
      auto res_lzs = lz_helper.find(i, 0xff, 2);
      auto res_l3 = lz_helper.find(i, 0x3fff, 3);
      dp0.update_lz(i, 2, 2, res_lzs, Constant<10>(), {lz2, 0, 0}, cost1);
      dp0.update_lz(i, 3, 3, res_lzs, Constant<11>(), {lz3, 0, 0}, cost1);
      dp0.update_lz(i, 3, 3, res_l3,  Constant<17>(), {lz3, 1, 0}, cost1);
      dp0.update_lz_matrix(i, ofs_tab, len_tab,
        [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
        [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
        0, cost1
      );
    }
    lz_helper.add_element(i);
  }

  if (dp1.optimal_cost() == dp1.infinite_cost) {
    // For example, this happens when the input is a de Bruijn sequence.
    throw std::runtime_error("This algorithm cannot compress the given data.");
  }

  auto commands = [&]{
    using command_type = sssp_solver<tag>::vertex_type;
    std::vector<command_type> commands;
    size_t curr = 1;
    ptrdiff_t adr = input.size();
    while (adr > 0 || (adr == 0 && curr == 1)) {
      command_type cmd;
      if (curr == 0) cmd = dp0[adr], curr = 1, assert(cmd.len > 0);
      else cmd = dp1[adr], curr = 0;
      adr -= cmd.len;
      commands.emplace_back(cmd);
    }
    assert(adr == 0 && curr == 0);
    std::ranges::reverse(commands);
    return commands;
  };

  using namespace data_type;
  writer_b8_l ret(4); ret.write<b1>(false);
  size_t adr = 0;
  for (auto&& cmd : commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type.tag) {
      case none: {
        ret.write<bnh>({3, 0});
      } break;

      case uncomp: {
        const size_t li = cmd.type.li;
        if (li == 0) ret.write<bnh>({3, cmd.len});
        else if (li == 1) ret.write<bnh, bnh>({4, 0b1110}, {4, cmd.len - 7});
        else ret.write<bnh, bnh>({4, 0b1111}, {10, cmd.len});
        ret.write<b8hn>({cmd.len, &input[adr]});
      } break;

      case lz2: {
        ret.write<bnh, bnh>({2, 0}, {8, d});
      } break;

      case lz3: {
        ret.write<bnh>({2, 1});
        if (d < 0x100) ret.write<bnh, bnh>({1, 0b1}, {8, d});
        else ret.write<bnh, bnh>({1, 0b0}, {14, d});
      } break;

      case lz: {
        const auto l = len_tab[cmd.type.li];
        const auto o = ofs_tab[cmd.type.oi];
        ret.write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
        ret.write<bnh>({o.bitlen, o.val + (d - o.min)});
      } break;

      default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp1.optimal_cost() + 33 == ret.bit_length());
  write32(ret.out, 0, input.size());

  const size_t s = std::min<size_t>(8, ret.size());
  uint64_t v = 1;
  for (size_t i = s - 1; i >= 4; --i) v = v << 8 | ret.out[i];
  v >>= 1;
  for (size_t i = 4; i < s; ++i) ret.out[i] = v & 0xff, v >>= 8;

  std::ranges::reverse(ret.out);
  return ret.out;
}

} // namespace sfc_comp
