#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> vortex_comp(std::span<const uint8_t> in) {
  enum CompType {
    none,
    uncomp_s, uncomp_m, uncomp_l,
    lz2, lz3, lz4, lz5, lz7, lzl
  };
  std::vector<uint8_t> input(in.begin(), in.end());
  std::reverse(input.begin(), input.end());

  lz_helper lz_helper(input);
  uncomp_helper u_helper(input.size(), 8);
  sssp_solver<CompType> dp0(input.size()), dp1(input.size());

  dp1[0].cost = dp1.infinite_cost;

  for (size_t i = 0; i <= input.size(); ++i) {
    const auto cost0 = dp0[i].cost;
    if (cost0 < dp0.infinite_cost) {
      dp1.update(i, 0, 0, Constant<3>(), none, cost0);
    }

    if (i == input.size()) break;

    if (cost0 < dp0.infinite_cost) {
      u_helper.update(i, cost0);
    }
    auto u0 = u_helper.find(i + 1, 1, 6);
    dp1.update_u(i + 1, u0.len, uncomp_s, u0.cost + 3);
    auto u1 = u_helper.find(i + 1, 7, 0x16);
    dp1.update_u(i + 1, u1.len, uncomp_m, u1.cost + 8);
    auto u2 = u_helper.find(i + 1, 0x17, 0x0fff);
    dp1.update_u(i + 1, u2.len, uncomp_l, u2.cost + 14);

    const auto cost1 = dp1[i].cost;
    if (cost1 < dp1.infinite_cost) {
      auto res_lzl = lz_helper.find_best(i, 0xffff);
      auto res_lzm = lz_helper.find_best(i, 0xfff);
      auto res_lzs = lz_helper.find_best(i, 0xff);
      auto res_l3 = lz_helper.find_best(i, 0x3fff);

      dp0.update_lz(i, 2, 2, res_lzs, Constant<10>(), lz2, cost1);
      dp0.update_lz(i, 3, 3, res_lzs, Constant<11>(), lz3, cost1);
      dp0.update_lz(i, 3, 3, res_l3,  Constant<17>(), lz3, cost1);

      static constexpr std::tuple<size_t, size_t, size_t, CompType> tab[] = {
        {4, 4, 12, lz4},
        {5, 6, 14, lz5},
        {7, 10, 16, lz7},
        {11, 0xff, 22, lzl}
      };

      for (const auto& tp : tab) {
        size_t beg, end, co; CompType ty;
        std::tie(beg, end, co, ty) = tp;
        dp0.update_lz(i, beg, end, res_lzs, Constant<0>(), ty, cost1 + co);
        dp0.update_lz(i, beg, end, res_lzm, Constant<4>(), ty, cost1 + co);
        dp0.update_lz(i, beg, end, res_lzl, Constant<7>(), ty, cost1 + co);
      }
    }

    lz_helper.add_element(i);
  }

  if (dp1.total_cost() == dp1.infinite_cost) {
    // For example, this happens when the input is a de Bruijn sequence.
    throw std::runtime_error("This algorithm cannot compress the given data.");
  }

  using command_type = sssp_solver<CompType>::vertex_type;

  auto commands = [&]{
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
    std::reverse(commands.begin(), commands.end());
    return commands;
  };

  using namespace data_type;
  writer_b ret; ret.write<b8ln_h, b1l>({32, 0}, 0);
  size_t adr = 0;
  for (auto&& cmd : commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
      case none: {
        ret.write<b8ln_h>({3, 0});
      } break;

      case uncomp_s: case uncomp_m: case uncomp_l: {
        if (cmd.type == uncomp_s) {
          ret.write<b8ln_h>({3, cmd.len});
        } else if (cmd.type == uncomp_m) {
          ret.write<b8ln_h, b8ln_h>({4, 0x0e}, {4, cmd.len - 7});
        } else {
          ret.write<b8ln_h, b8ln_h>({4, 0x0f}, {10, cmd.len});
        }
        for (size_t i = 0; i < cmd.len; ++i) {
          ret.write<b8ln_h>({8, input[adr + i]});
        }
      } break;

      case lz2: case lz3:
      case lz4: case lz5: case lz7: case lzl: {
        if (cmd.type == lz2) {
          ret.write<b8ln_h>({2, 0});
          ret.write<b8ln_h>({8, d});
        } else if (cmd.type == lz3) {
          ret.write<b8ln_h>({2, 1});
          if (d < 0x100) {
            ret.write<b8ln_h>({1, 1});
            ret.write<b8ln_h>({8, d});
          } else {
            ret.write<b8ln_h>({1, 0});
            ret.write<b8ln_h>({14, d});
          }
        } else {
          if (cmd.type == lz4) {
            ret.write<b8ln_h>({2, 2});
          } else {
            ret.write<b8ln_h>({2, 3});
            if (cmd.type == lz5) {
              ret.write<b8ln_h>({2, cmd.len - 5});
            } else if (cmd.type == lz7) {
              ret.write<b8ln_h>({2, 2});
              ret.write<b8ln_h>({2, cmd.len - 7});
            } else {
              ret.write<b8ln_h>({2, 3});
              ret.write<b8ln_h>({8, cmd.len});
            }
          }
          if (d < 0x100) {
            ret.write<b8ln_h>({2, 3});
            ret.write<b8ln_h>({8, d});
          } else if (d < 0x1000) {
            ret.write<b8ln_h>({2, 2});
            ret.write<b8ln_h>({12, d});
          } else {
            ret.write<b8ln_h>({1, 0});
            ret.write<b8ln_h>({16, d});
          }
        }
      }
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp1.total_cost() + 33 == 8 * ret.size() - ret.bit);
  write32(ret.out, 0, input.size());
  uint64_t v = 1;
  for (size_t i = std::min<size_t>(8, ret.out.size()) - 1; i >= 4; --i) {
    v = v * 256 + ret.out[i];
  }
  v >>= 1;
  for (size_t i = 4; v > 0; ++i, v >>= 8) {
    ret.out[i] = v;
  }
  std::reverse(ret.out.begin(), ret.out.end());
  return ret.out;
}

} // namespace sfc_comp
