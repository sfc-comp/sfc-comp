#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> mahoujin_guru_guru_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum method { uncomp, uncomp0, prev2, rle0 };
  using tag = tag_l<method>;

  static constexpr auto rle_lens = create_array<vrange, 8>([](size_t k) {
    return vrange((1 << k) + 1, std::min(255, 2 << k), 2 * k + 1, 0);
  });

  sssp_solver<tag> dp(input.size());

  size_t rlen = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    rlen = encode::run_length(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, 1, 1, Constant<2>(), {uncomp0, 0});
      const auto cost = dp[i].cost;
      for (size_t k = 0; k < rle_lens.size(); ++k) {
        dp.update(i, rle_lens[k].min, rle_lens[k].max, rlen,
                  Constant<3>(), {rle0, k}, cost + rle_lens[k].bitlen);
      }
    } else {
      dp.update(i, 1, 1, Constant<9>(), {uncomp, 0});
    }
    if (i >= 2 && input[i] == input[i - 2]) dp.update(i, 1, 1, Constant<3>(), {prev2, 0});
  }

  using namespace data_type;
  writer_b8_h ret(10);
  writer raw;

  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1>(false);
      raw.write<d8>(input[adr]);
    } break;
    case uncomp0: {
      ret.write<b1, b1>(true, false);
    } break;
    case rle0: {
      const auto k = cmd.type.li;
      ret.write<b1, b1, b1>(true, true, false);
      ret.write<bnh>({k, (size_t(1) << k) - 1});
      ret.write<b1>(false);
      const size_t l = cmd.len - (rle_lens[0].min - 1);
      ret.write<bnh>({k, l ^ (1 << k)});
    } break;
    case prev2: {
      ret.write<b1, b1, b1>(true, true, true);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.total_cost() + 10 * 8 == raw.size() * 8 + ret.bit_length());

  write32b(ret.out, 0, 0x4e4d5030); // NMP0
  write16(ret.out, 4, raw.size());
  write16(ret.out, 6, ret.size() - 10);
  write16(ret.out, 8, input.size());
  std::copy(raw.out.begin(), raw.out.end(), std::back_inserter(ret.out));

  return ret.out;
}

} // namespace sfc_comp
