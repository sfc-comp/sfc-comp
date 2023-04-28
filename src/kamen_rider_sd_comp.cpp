#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> kamen_rider_sd_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);

  enum CompType { uncomp, lzs, lz, rle0, rle1 };

  static constexpr size_t rle_max_lens[2] = {9, 2};
  static constexpr size_t tab_max_size = 0x11;

  struct thres_size {
    bool operator == (const thres_size& rhs) const {
      return thres == rhs.thres && size == rhs.size;
    }
    size_t thres;
    size_t size;
  };

  struct pre_rle : public std::vector<uint8_t> {
    pre_rle() : thres(0) {}
    pre_rle(const std::vector<uint8_t>& in, size_t thres) : std::vector<uint8_t>(in), thres(thres) {}
    bool oversized() const { return thres > 2 || size() > tab_max_size; }
    bool valid() const { return !oversized() && ((thres == 2 && size() >= 2) || (thres == size())); }
    size_t thres;
  };

  const auto pre_20 = pre_rle({
    0x00, 0xff,
    0x80, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc, 0xfe, 0x7f, 0x3f, 0x1f, 0x0f, 0x07, 0x03, 0x01, 0x18
  }, 2);
  const auto pre_40 = pre_rle({
    0x9e, 0xc9,
    0x92, 0xcb, 0xb9, 0xb0, 0xbc, 0xaf, 0xbb, 0xa8, 0xa6, 0xa4, 0xb5, 0xad, 0xb7, 0xca, 0xa9
  }, 2);

  lz_helper lz_helper(input);
  std::vector<encode::lz_data> lzs_memo(input.size()), lz_memo(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    lzs_memo[i] = lz_helper.find_best(i, 0x20);
    lz_memo[i] = lz_helper.find_best(i, 0x1000);
    lz_helper.add_element(i);
  }

  auto pre = [&]{
    std::array<size_t, 256> gain = {};
    const size_t max_len = rle_max_lens[0];
    for (size_t i = 0; i < input.size(); ) {
      size_t l = std::min<size_t>(max_len, encode::run_length(input, i, 0));
      if (l >= 2) gain[input[i]] += l - 1;
      i += l;
    }
    std::vector<uint8_t> vals;
    for (size_t i = 0; i < gain.size(); ++i) if (gain[i] > 0) vals.push_back(i);
    std::sort(vals.begin(), vals.end(), [gain](size_t a, size_t b) { return gain[a] > gain[b]; });
    return pre_rle(vals, vals.size());
  }();

  struct result {
    std::vector<uint8_t> compressed;
    pre_rle pre;
  };

  const auto compress = [&](const pre_rle& pre, const size_t header_val, bool output, const thres_size limit) -> result {
    std::array<std::array<ptrdiff_t, 256>, 2> ind;
    for (auto& v : ind) std::fill(v.begin(), v.end(), - 1);
    assert(pre.thres <= pre.size());
    for (size_t k = 0; k < pre.thres; ++k) ind[0][pre[k]] = k;
    for (size_t k = pre.thres; k < pre.size(); ++k) ind[1][pre[k]] = k - pre.thres;

    sssp_solver<CompType> dp(input.size());

    size_t rlen = 0;
    for (size_t i = 0; i < input.size(); ++i) {
      dp.update(i, 1, 1, Constant<9>(), uncomp);
      dp.update_lz(i, 2, 6, lzs_memo[i], Constant<9>(), lzs);
      dp.update_lz(i, 3, 6, lz_memo[i], Constant<17>(), lz);
      rlen = encode::run_length(input, i, rlen);
      if (rlen <= 1) continue;
      for (size_t l = 0; l < 2; ++l) {
        if (auto k = ind[l][input[i]]; k >= 0) {
          dp.update_lz(i, 2, rle_max_lens[l], {k, rlen}, Constant<9>(), (l == 0) ? rle0 : rle1);
        }
      }
    }

    using namespace data_type;
    const auto commands = dp.commands();

    if (!output) {
      std::array<size_t, 256> gain = {};
      for (const auto& cmd : commands) {
        if (cmd.type == rle0) gain[pre[cmd.val()]] += cmd.len - 1;
        else if (cmd.type == rle1) gain[pre[pre.thres + cmd.val()]] += cmd.len - 1;
      }
      const auto pick = [](const std::span<const uint8_t> vals, const std::span<const size_t> gain, size_t nsize) {
        std::vector<uint8_t> ret;
        for (auto v : vals) if (gain[v] > 0) ret.push_back(v);
        nsize = std::min(nsize, ret.size());
        std::partial_sort(ret.begin(), ret.begin() + nsize, ret.end(),
                          [&](size_t a, size_t b) { return gain[a] > gain[b]; });
        ret.resize(nsize);
        // std::sort(ret.begin(), ret.end());
        return ret;
      };
      assert(limit.size >= limit.thres);
      const auto first = pick(pre, gain, limit.thres);
      for (auto v : first) gain[v] = 0;
      auto ret = pre_rle(first, first.size());

      const auto second = pick(pre, gain, limit.size - first.size());
      std::copy(second.begin(), second.end(), std::back_inserter(ret));
      return {{}, ret};
    } else {
      if (!pre.valid()) throw std::logic_error("Broken RLE table.");
      writer_b ret;
      if (header_val >= 0x20) {
        ret.write<d8>(header_val);
      } else {
        ret.write<d8, d8n>(pre.size(), {pre.size(), pre.data()});
      }
      size_t adr = 0;
      for (const auto& cmd : commands) {
        const size_t d = adr - cmd.lz_ofs;
        switch (cmd.type) {
        case uncomp: ret.write<b1h, d8>(false, input[adr]); break;
        case lzs: ret.write<b1h, d8>(true, 0x00 | (cmd.len - 2) << 5 | (d - 1)); break;
        case rle0: ret.write<b1h, d8>(true, 0xa0 | (cmd.val() << 3) | (cmd.len - 2)); break;
        case rle1: ret.write<b1h, d8>(true, 0xb0 | cmd.val()); break;
        case lz: ret.write<b1h, d16b>(true, 0xc000 | (cmd.len - 3) << 12 | (d - 1)); break;
        default: break;
        }
        adr += cmd.len;
      }
      ret.write<b1h, d8>(true, 0xbf);
      assert(adr == input.size());
      if (header_val >= 0x20) {
        assert(dp.total_cost() + 1 + (1 + 1) * 8 == ret.bit_length());
      } else {
        assert(dp.total_cost() + 1 + (1 + pre.size() + 1) * 8 == ret.bit_length());
      }
      return {ret.out, {}};
    }
  };

  const std::vector<size_t> coeffs = {128, 64, 32, 16, 8, 0};
  std::vector<uint8_t> best;

  thres_size last_limit = {pre.thres, pre.size()};
  for (size_t lv = 0; lv < coeffs.size(); ++lv) {
    const size_t c = coeffs[lv];
    const size_t nsize = (pre.size() <= tab_max_size) ? pre.size()
                                                      : tab_max_size + (pre.size() - tab_max_size) * c / 256;
    const size_t nthres = (nsize <= 2) ? nsize : (2 + (nsize - 2) / tab_max_size);
    const auto limit = thres_size({nthres, nsize});
    if (limit == last_limit) continue;
    last_limit = limit;
    pre = compress(pre, 0, false, limit).pre;
  }
  pre = compress(pre, 0, false, last_limit).pre; // remove unused rles.
  best = std::move(compress(pre, 0, true, {}).compressed);
  assert(best.size() > 0);

  if (auto res = compress(pre_20, 0x20, true, {}).compressed; res.size() < best.size()) {
    best = std::move(res);
  }
  if (auto res = compress(pre_40, 0x40, true, {}).compressed; res.size() < best.size()) {
    best = std::move(res);
  }
  return best;
}

} // namespace sfc_comp
