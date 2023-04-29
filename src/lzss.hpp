#pragma once

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

template <class Writer, typename InitFunc, typename LzEncoding>
std::vector<uint8_t> lzss(
    std::span<const uint8_t> in,
    const size_t pad, InitFunc init,
    const size_t lz_max_ofs,
    const size_t lz_min_len, const size_t lz_max_len,
    const size_t header_size,
    const bool uncomp_flag, LzEncoding lz_enc) {

  enum CompType { uncomp, lz };
  std::vector<uint8_t> input(in.size() + pad);
  std::copy(in.begin(), in.end(), input.begin() + pad);

  init(input);

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);
  for (size_t i = 0; i < pad; ++i) dp[i + 1].cost = 0;

  for (size_t i = pad; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, lz_max_ofs);
    dp.update_lz(i, lz_min_len, lz_max_len, res_lz, Constant<17>(), lz);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  Writer ret(header_size);

  size_t adr = pad;
  for (const auto cmd : dp.commands(pad)) {
    switch (cmd.type) {
    case uncomp: ret.template write<b1, d8>(uncomp_flag, input[adr]); break;
    case lz: ret.template write<b1, d16>(!uncomp_flag, lz_enc(adr, cmd.lz_ofs, cmd.len)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.total_cost() + header_size * 8 == ret.bit_length());
  assert(adr - pad == in.size());
  return ret.out;
}

} // namespace sfc_comp
