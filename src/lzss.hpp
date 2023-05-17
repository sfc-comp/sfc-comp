#pragma once

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

template <class Writer, typename InitFunc, typename LzEncoding>
requires std::derived_from<Writer, writer> &&
         std::invocable<InitFunc, std::span<uint8_t>> &&
         std::invocable<LzEncoding, size_t, size_t, size_t>
std::vector<uint8_t> lzss(
    std::span<const uint8_t> in,
    const size_t pad, InitFunc&& init,
    const size_t lz_max_ofs, const size_t lz_min_len, const size_t lz_max_len,
    const size_t header_size,
    const bool uncomp_b, LzEncoding&& lz_enc) {

  enum tag { uncomp, lz };
  std::vector<uint8_t> input(in.size() + pad);
  std::ranges::copy(in, input.begin() + pad);

  init(input);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size()); auto c0 = dp.template c<0>(lz_max_len);

  for (size_t i = input.size(); i-- > pad; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, uncomp);
    auto res_lz = lz_helper.find(i, lz_max_ofs, lz_min_len);
    dp.update(i, lz_min_len, lz_max_len, res_lz, c0, 17, lz);
    c0.update(i);
  }

  using namespace data_type;
  Writer ret(header_size);

  size_t adr = pad;
  for (const auto& cmd : dp.optimal_path(pad)) {
    switch (cmd.type) {
    case uncomp: ret.template write<b1, d8>(uncomp_b, input[adr]); break;
    case lz: ret.template write<b1, d16>(!uncomp_b, lz_enc(adr, cmd.lz_ofs(), cmd.len)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.optimal_cost(pad) + header_size * 8 == ret.bit_length());
  assert(adr == input.size());
  return ret.out;
}

} // namespace sfc_comp
