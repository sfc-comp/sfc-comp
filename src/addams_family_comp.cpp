#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> addams_family_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x800000);

  enum CompType {
    uncomp, lzs, lzm, lzl
  };
  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x400);
    dp.update_lz(i, 2, 2, res_lz, Constant<12>(), lzs);
    dp.update_lz(i, 3, 9, res_lz, Constant<15>(), lzm);
    dp.update_lz(i, 10, 264, res_lz, Constant<23>(), lzl);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret;
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type) {
    case uncomp: {
      ret.write<b1, bnh>(false, {8, input[adr]});
    } break;
    case lzs:
    case lzm:
    case lzl: {
      ret.write<b1>(true);
      if (cmd.type == lzs) ret.write<bnh>({1, 0});
      else if (cmd.type == lzm) ret.write<bnh>({4, 8 + (cmd.len - 2)});
      else ret.write<bnh, bnh>({4, 8}, {8, cmd.len - 9});
      ret.write<bnh>({10, adr - cmd.lz_ofs - 1});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh>({13, 0x1800});
  assert(adr == input.size());
  assert(dp.total_cost() + 13 == ret.bit_length());
  return ret.out;
}

std::vector<uint8_t> jurassic_park_comp(std::span<const uint8_t> input) {
  auto ret = addams_family_comp(input);
  if (ret.size() & 1) ret.push_back(0);
  for (size_t i = 0; i < ret.size(); i += 2) std::swap(ret[i], ret[i + 1]);
  return ret;
}

} // namespace sfc_comp
