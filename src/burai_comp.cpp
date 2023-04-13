#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> burai_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum CompType {
    uncomp,
    lzs1, lzs2, lzs3,
    lzl2, lzl3, lzl4, lzl
  };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  struct codeword {
    size_t bits;
    size_t val;
  };
  std::array<codeword, 256> codewords = {};
  for (size_t i = 0; i < 256; ++i) codewords[i] = {10, 0x300 | i};
  for (size_t i = 0; i < 16; ++i) codewords[i] = {7, 0x40 | i};
  codewords[0x10] = {5, 0x14}; codewords[0x30] = {5, 0x15};
  codewords[0x80] = {5, 0x16}; codewords[0xff] = {5, 0x17};

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, [&](size_t) { return codewords[input[i]].bits; }, uncomp);
    auto res_lz4 = lz_helper.find_best(i, 0x10);
    auto res_lz8 = lz_helper.find_best(i, 0x100);

    dp.update_lz(i, 1, 1, res_lz4, Constant<6>(), lzs1);
    dp.update_lz(i, 2, 2, res_lz4, Constant<7>(), lzs2);
    dp.update_lz(i, 3, 3, res_lz4, Constant<8>(), lzs3);
    dp.update_lz(i, 2, 2, res_lz8, Constant<13>(), lzl2);
    dp.update_lz(i, 3, 3, res_lz8, Constant<14>(), lzl3);
    dp.update_lz(i, 4, 4, res_lz8, Constant<15>(), lzl4);
    dp.update_lz(i, 5, 20, res_lz8, Constant<20>(), lzl);

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const auto d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: {
      const auto c = codewords[input[adr]];
      ret.write<b8hn_h>({c.bits, c.val});
    } break;
    case lzs1:
    case lzs2:
    case lzs3: {
      if (cmd.type == lzs1) ret.write<b8hn_h>({2, 0});
      else if (cmd.type == lzs2) ret.write<b8hn_h>({3, 2});
      else ret.write<b8hn_h>({4, 6});
      ret.write<b8hn_h>({4, d});
    } break;
    case lzl2:
    case lzl3:
    case lzl4:
    case lzl: {
      if (cmd.type == lzl2) ret.write<b8hn_h>({5, 0x0e});
      else if (cmd.type == lzl3) ret.write<b8hn_h>({6, 0x1e});
      else if (cmd.type == lzl4) ret.write<b8hn_h>({7, 0x3e});
      else ret.write<b8hn_h, b8hn_h>({8, 0x7e}, {4, cmd.len - 5});
      ret.write<b8hn_h>({8, d});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<b8hn_h>({8, 0x7f});
  assert((dp.total_cost() + 8 + 7) / 8 == ret.size());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
