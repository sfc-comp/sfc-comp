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
  std::copy(in.begin(), in.end(), input.begin() + pad);

  sssp_solver<tag> dp(input.size(), pad);

  size_t rleni = encode::run_length(input, pad - 1, 0, 1);
  std::array<size_t, dists.size()> res_lz = {};
  for (size_t i = pad; i < input.size(); ++i) {
    dp.update(i, 1, 0x20, Linear<1, 1>(), uncomp);

    if (rleni >= 2) dp.update(i, 1, 0x20, rleni - 1, Constant<1>(), inc);
    rleni = encode::run_length(input, i, rleni, 1);

    for (size_t k = 0; k < dists.size(); ++k) {
      res_lz[k] = encode::lz_dist(input, i, dists[k], res_lz[k]);
      dp.update_lz(i, 1, 0x20, {k + 2, res_lz[k]}, Constant<1>(), lz);
    }
  }

  using namespace data_type;
  writer ret;
  size_t adr = pad;
  for (const auto cmd : dp.commands(pad)) {
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>((0 << 5) | (cmd.len - 1), {cmd.len, &input[adr]}); break;
    case inc: ret.write<d8>((1 << 5) | (cmd.len - 1)); break;
    case lz: ret.write<d8>((cmd.val() << 5) | (cmd.len - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(dp.optimal_cost() == ret.size());
  assert(adr == input.size());

  return ret.out;
}

} // namespace sfc_comp
