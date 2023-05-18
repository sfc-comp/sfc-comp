#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> knights_of_the_round_comp(std::span<const uint8_t> in) {
  check_size(in.size(), 0, 0xffff);

  enum tag { uncomp, uncompl, lzs, lz, pre0, pre1, pre2, pre3, pre };

  static constexpr size_t pad = 0x12;
  std::vector<uint8_t> input(in.size() + pad, 0x20);
  std::ranges::copy(in, input.begin() + pad);

  lz_helper lz_helper(input);
  std::vector<std::array<encode::lz_data, 2>> lz_memo(input.size());
  for (size_t i = 0; i < pad; ++i) lz_helper.add_element(i);
  for (size_t i = pad; i < input.size(); ++i) {
    lz_memo[i][0] = lz_helper.find(i, 0x100, 2);
    lz_memo[i][1] = lz_helper.find(i, 0x1000, 3);
    lz_helper.add_element(i);
  }

  static constexpr size_t iter_total = 5;
  static constexpr auto cand_size = std::to_array<size_t, iter_total>({32, 16, 12, 8, 4});
  static constexpr size_t dummy_bits = 6;
  static_assert(iter_total == 5 && cand_size.back() == 4);

  const auto freq = utility::freq_u8(in);
  std::vector<uint8_t> cands;
  for (const auto v : utility::k_most<uint8_t, cand_size[0]>(freq)) cands.push_back(v);

  for (size_t iter = 0; iter < iter_total; ++iter) {
    solver<tag> dp(input.size());
    auto c0 = dp.c<0>(3 + 0x0f);
    auto c8 = dp.c<8>(0x10 + 0xff);

    std::vector<ptrdiff_t> ind(256, -1);
    for (size_t k = 0; k < cands.size(); ++k) ind[cands[k]] = k;

    for (size_t i = input.size(); i-- > pad; ) {
      dp.update(i, 1, 1, c8, 1, uncomp);
      dp.update(i, 0x10, 0x10 + 0xff, c8, 15, uncompl);

      dp.update(i, 2, 2 + 0x07, lz_memo[i][0], c0, 15, lzs);
      dp.update(i, 3, 3 + 0x0f, lz_memo[i][1], c0, 18, lz);

      const auto v = input[i];
      if (auto k = ind[v]; k >= 0) {
        static constexpr tag tags[4] = {pre0, pre1, pre2, pre3};
        static constexpr size_t bits[4] = {3, 5, 6, 7};
        if (size_t(k) < iter) dp.update(i, 1, bits[k], tags[k]);
        else dp.update(i, 1, dummy_bits, pre);
      }
      c0.update(i); c8.update(i);
    }

    if (iter + 1 < iter_total) {
      std::array<size_t, 256> counter = {};
      size_t adr = pad;
      for (const auto& cmd : dp.optimal_path(adr)) {
        if (cmd.type == pre) counter[input[adr]] += 1;
        adr += cmd.len;
      }
      assert(adr == input.size());
      const size_t nsize = cand_size[iter + 1];
      std::partial_sort(cands.begin() + iter, cands.begin() + nsize, cands.end(),
        [&](size_t a, size_t b) { return counter[a] > counter[b]; });
      cands.resize(nsize);
    } else {
      using namespace data_type;
      writer_b8_l ret(2);
      writer raw(6);

      size_t adr = pad;
      for (const auto& cmd : dp.optimal_path(pad)) {
        const size_t d = adr - cmd.lz_ofs();
        switch (cmd.type) {
        case uncomp: {
          ret.write<b1>(true);
          raw.write<d8>(input[adr]);
        } break;
        case lz: {
          size_t v = (cmd.lz_ofs() - 0x12 - pad) & 0xfff;
          ret.write<bnh>({2, 1});
          raw.write<d16>((v & 0x0f00) << 4 | (cmd.len - 3) << 8 | (v & 0x00ff));
        } break;
        case pre0: ret.write<bnh>({3, 1}); break;
        case pre1: ret.write<bnh>({5, 1}); break;
        case pre2: ret.write<bnh>({6, 1}); break;
        case pre3: ret.write<bnh>({7, 1}); break;
        case lzs: {
          ret.write<bnh>({7, 8 + (cmd.len - 2)});
          raw.write<d8>(0x100 - d);
        } break;
        case uncompl: {
          ret.write<bnh>({7, 0});
          raw.write<d8, d8n>(cmd.len - 0x10, {cmd.len, &input[adr]});
        } break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      write16(ret.out, 0, ret.size() - 2);
      write16(raw.out, 0, in.size());
      std::copy_n(cands.begin(), 4, raw.out.begin() + 2);
      assert(adr == input.size());
      assert(dp.optimal_cost(pad) + 8 * 8 == raw.size() * 8 + ret.bit_length());
      ret.extend(raw);
      return ret.out;
    }
  }
  throw std::logic_error("iter_total == 0");
}

} // namespace sfc_comp
