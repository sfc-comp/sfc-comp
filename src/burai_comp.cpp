#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> burai_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum tag {
    uncomp,
    lzs1, lzs2, lzs3,
    lzl2, lzl3, lzl4, lzl
  };

  static constexpr auto codewords = create_array<encode::codeword, 256>([](size_t i) -> encode::codeword {
    if (i < 0x10) return {7, 0x40 | i};
    if (i == 0x10) return {5, 0x14};
    if (i == 0x30) return {5, 0x15};
    if (i == 0x80) return {5, 0x16};
    if (i == 0xff) return {5, 0x17};
    return {10, 0x300 | i};
  });

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, [&](size_t) { return codewords[input[i]].bitlen; }, uncomp);
    auto res_lz4 = lz_helper.find(i, 0x10, 1);
    dp.update_lz(i, 1, 1, res_lz4, Constant<6>(), lzs1);
    dp.update_lz(i, 2, 2, res_lz4, Constant<7>(), lzs2);
    dp.update_lz(i, 3, 3, res_lz4, Constant<8>(), lzs3);

    auto res_lz8 = lz_helper.find(i, 0x100, 2);
    dp.update_lz(i, 2, 2, res_lz8, Constant<13>(), lzl2);
    dp.update_lz(i, 3, 3, res_lz8, Constant<14>(), lzl3);
    dp.update_lz(i, 4, 4, res_lz8, Constant<15>(), lzl4);
    dp.update_lz(i, 5, 20, res_lz8, Constant<20>(), lzl);

    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    const auto d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: {
      ret.write<bnh>(codewords[input[adr]]);
    } break;
    case lzs1:
    case lzs2:
    case lzs3: {
      if (cmd.type == lzs1) ret.write<bnh>({2, 0});
      else if (cmd.type == lzs2) ret.write<bnh>({3, 2});
      else ret.write<bnh>({4, 6});
      ret.write<bnh>({4, d});
    } break;
    case lzl2:
    case lzl3:
    case lzl4:
    case lzl: {
      if (cmd.type == lzl2) ret.write<bnh>({5, 0x0e});
      else if (cmd.type == lzl3) ret.write<bnh>({6, 0x1e});
      else if (cmd.type == lzl4) ret.write<bnh>({7, 0x3e});
      else ret.write<bnh, bnh>({8, 0x7e}, {4, cmd.len - 5});
      ret.write<bnh>({8, d});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({8, 0x7f});
  assert(adr == input.size());
  assert(dp.optimal_cost() + 8 == ret.bit_length());
  return ret.out;
}

} // namespace sfc_comp
