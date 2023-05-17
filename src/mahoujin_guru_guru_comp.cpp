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

  solver<tag> dp(input.size());
  auto c0 = dp.c<0>(rle_lens.back().max);

  size_t rlen = 0;
  for (size_t i = input.size(); i-- > 0; ) {
    rlen = encode::run_length_r(input, i, rlen);
    if (input[i] == 0) {
      dp.update(i, 1, 2, {uncomp0, 0});
      dp.update(i, rle_lens, rlen, c0, 3, [&](size_t li) -> tag { return {rle0, li}; });
    } else {
      dp.update(i, 1, 9, {uncomp, 0});
    }
    if (i >= 2 && input[i] == input[i - 2]) dp.update(i, 1, 3, {prev2, 0});
    c0.update(i);
  }

  using namespace data_type;
  writer_b8_h ret(10);
  writer raw;

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      ret.write<b1>(false);
      raw.write<d8>(input[adr]);
    } break;
    case uncomp0: {
      ret.write<b1, b1>(true, false);
    } break;
    case rle0: {
      ret.write<b1, b1, b1>(true, true, false);
      ret.write<bnh>({li, low_bits_mask(li)});
      ret.write<b1>(false);
      const size_t l = cmd.len - (rle_lens[0].min - 1);
      ret.write<bnh>({li, l ^ (1 << li)});
    } break;
    case prev2: {
      ret.write<b1, b1, b1>(true, true, true);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  assert(adr == input.size());
  assert(dp.optimal_cost() + 10 * 8 == raw.size() * 8 + ret.bit_length());

  write32b(ret.out, 0, 0x4e4d5030); // NMP0
  write16(ret.out, 4, raw.size());
  write16(ret.out, 6, ret.size() - 10);
  write16(ret.out, 8, input.size());
  ret.extend(raw);

  return ret.out;
}

} // namespace sfc_comp
