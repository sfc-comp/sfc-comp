#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> super_bomberman_5_comp(std::span<const uint8_t> in) {
  // This compression algorithm works for every input,
  // but the decompressor only accepts inputs whose size is divisible by 2.
  check_size(in.size(), 2, 0x8100);

  enum tag { uncomp, inc, lz };

  static constexpr auto dists = std::to_array<size_t>({1, 2, 4, 16, 32, 64});
  static constexpr size_t pad = 0x100;

  std::vector<uint8_t> input(in.size() + pad);
  std::ranges::copy(in, input.begin() + pad);

  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(0x20);
  auto c1 = dp.c<1>(0x20);

  size_t rleni = 0;
  std::array<size_t, dists.size()> res_lz = {};

  for (size_t i = input.size(); i-- > pad; ) {
    dp.update(i, 1, 0x20, c1, 1, uncomp);

    rleni = encode::lz_dist_r(input, i, 1, rleni, 1);
    dp.update(i, 1, 0x20, rleni, c0, 1, inc);

    for (size_t k = 0; k < dists.size(); ++k) {
      res_lz[k] = encode::lz_dist_r(input, i, dists[k], res_lz[k]);
    }
    const size_t best_k = std::max_element(res_lz.begin(), res_lz.end()) - res_lz.begin();
    dp.update(i, 1, 0x20, {best_k + 2, res_lz[best_k]}, c0, 1, lz);

    c0.update(i); c1.update(i);
  }

  using namespace data_type;
  writer ret;
  size_t adr = pad;
  for (const auto& cmd : dp.optimal_path(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>((0 << 5) | (cmd.len - 1), {cmd.len, &input[adr]}); break;
    case inc: ret.write<d8>((1 << 5) | (cmd.len - 1)); break;
    case lz: ret.write<d8>((cmd.arg << 5) | (cmd.len - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.optimal_cost(pad) == ret.size());
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
