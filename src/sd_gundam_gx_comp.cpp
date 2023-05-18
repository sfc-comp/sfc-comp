#include <numeric>

#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

#include "huffman.hpp"

namespace sfc_comp {

namespace {

struct shannon_fano {
  std::vector<size_t> words;
  std::vector<encode::codeword> codewords;
  std::vector<size_t> bits;
};

shannon_fano shannon_fano_encode(std::span<const size_t> counts, bool add_header_cost = false) {
  using offset_array = std::array<size_t, 8>;

  struct sf_data {
    sf_data(size_t prefix_len) : prefix_len(prefix_len), bits(1 << prefix_len) {}
    size_t cost = std::numeric_limits<size_t>::max();
    size_t prefix_len;
    std::vector<size_t> bits;
    offset_array offsets;
    size_t header_cost;
  };

  const auto words_encoding_cost = [&counts](std::span<const size_t> words) {
    static constexpr auto encode_costs = std::to_array<size_t>({
      2, 3, 3, 6, 6, 6, 6, 6,
      6, 6, 6, 6, 6, 6, 6, 6,
      6, 6, 6
    });
    size_t prev_word = counts.size();
    size_t ret = 0;
    for (size_t j = 0; j < words.size(); ++j) {
      size_t word = words[j];
      if (word <= prev_word || word >= prev_word + 0x14) ret += 11;
      else ret += encode_costs[word - prev_word - 1];
      prev_word = word;
    }
    return ret;
  };

  const auto bucket_sort = [&](std::span<const size_t> words, const offset_array& offsets,
      std::span<const size_t> bits, std::span<size_t> dest) {
    offset_array ofs; std::ranges::copy(offsets, ofs.begin());
    for (const auto w : words) dest[ofs[bits[w] - 1]++] = w;
  };

  std::vector<size_t> words, ordered_words;
  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] <= 0) continue;
    words.push_back(i);
    ordered_words.push_back(i);
  }
  if (words.size() == 0) return shannon_fano();

  std::sort(words.begin(), words.end(), [&](const size_t a, const size_t b) {
    return counts[a] > counts[b] || (counts[a] == counts[b] && a < b); }
  );

  std::vector<size_t> cumu(words.size());
  for (size_t i = 0; i < cumu.size(); ++i) cumu[i] = counts[words[i]];
  for (size_t i = cumu.size() - 1; i > 0; --i) cumu[i - 1] += cumu[i];

  offset_array curr_offsets = {};
  std::array<size_t, 8> curr_bits;
  std::vector<size_t> curr_words(words.size());
  std::vector<size_t> word_bits(counts.size());

  std::function<sf_data(size_t, size_t, size_t, size_t, size_t)> solve =
      [&](size_t now, size_t min_bit, size_t num_words, size_t bit_size, size_t cost) {

    const size_t size = 1 << min_bit;

    sf_data ret(min_bit);

    for (size_t b = bit_size; b <= 8; ++b) {
      const size_t count = 1 << b;
      curr_offsets[b - 1] = num_words;

      size_t least_num_words = num_words + count * (size - now);
      if (least_num_words >= words.size()) {
        if (cost >= ret.cost) break;
        size_t encode_cost = cost;
        if (add_header_cost) {
          for (size_t j = num_words; j < words.size(); ++j) word_bits[words[j]] = b;
          bucket_sort(ordered_words, curr_offsets, word_bits, curr_words);
          encode_cost += words_encoding_cost(curr_words);
        }
        if (encode_cost < ret.cost) {
          ret.cost = encode_cost;
          ret.header_cost = encode_cost - cost;
          std::fill_n(curr_bits.begin() + now, size - now, b);
          std::copy_n(curr_bits.begin(), size, ret.bits.begin());
          std::ranges::copy(curr_offsets, ret.offsets.begin());
        }
        break;
      }

      for (size_t i = now; i < size - 1; ++i) {
        const size_t curr = num_words + count * (i - now), next = curr + count;
        const size_t ncost = cost + cumu[next];
        if (ncost >= ret.cost) break;
        curr_bits[i] = b;
        if (add_header_cost) {
          for (size_t j = curr; j < next; ++j) word_bits[words[j]] = b;
        }
        auto res = solve(i + 1, min_bit, next, b + 1, ncost);
        if (res.cost < ret.cost) ret = std::move(res);
      }
      cost += cumu[num_words];
      if (cost >= ret.cost) break;
    }
    return ret;
  };

  sf_data best_4 = solve(0, 2, 0, 1, 3 * cumu[0] + 4 * 3);
  sf_data best_8 = solve(0, 3, 0, 1, 4 * cumu[0] + 8 * 3);
  const sf_data best = std::move((best_4.cost < best_8.cost) ? best_4 : best_8);

  size_t curr = 0;
  for (const auto b : best.bits) {
    size_t next = std::min(words.size(), curr + (1 << b));
    for (size_t j = curr; j < next; ++j) word_bits[words[j]] = b;
    curr = next;
  }
  bucket_sort(ordered_words, best.offsets, word_bits, curr_words);

  std::vector<encode::codeword> codewords(counts.size(), {-1, 0});
  curr = 0;
  for (size_t i = 0; i < best.bits.size(); ++i) {
    size_t b = best.bits[i];
    size_t next = std::min(words.size(), curr + (1 << b));
    ptrdiff_t bitlen = best.prefix_len + b;
    for (size_t j = curr; j < next; ++j) codewords[curr_words[j]] = {bitlen, (i << b) | (j - curr)};
    curr = next;
  }
  return shannon_fano { .words = std::move(curr_words),
                        .codewords = std::move(codewords),
                        .bits = std::move(best.bits) };
}

} // namespace

std::vector<uint8_t> sd_gundam_gx_comp_core(std::span<const uint8_t> input, const bool dsp3) {
  if (dsp3) {
    check_size(input.size(), 1, 0xffff);
  } else {
    check_size(input.size(), 0, 0xffff);
    if (input.size() == 0) {
      return std::vector<uint8_t>(2, 0);
    }
  }

  enum method { uncomp, lz };
  using tag = tag_o<method>;

  static constexpr size_t lz_min_len = 3;
  static constexpr size_t lz_max_len = lz_min_len + 0xff;
  static constexpr size_t sf_offset = 0x0100 - lz_min_len;
  static constexpr auto lz_lens = create_array<size_t, lz_max_len - lz_min_len + 1>([&](size_t i) {
    return lz_min_len + i;
  });

  const size_t ofs_bits_s = (dsp3) ? 8 : 9;
  const size_t ofs_bits_l = ofs_bits_s + 4;
  const auto ofs_tab = to_vranges({
    {size_t(1),                     ofs_bits_s + 1, size_t(0b0 << ofs_bits_s)},
    {size_t(1) + (1 << ofs_bits_s), ofs_bits_l + 1, size_t(0b1 << ofs_bits_l) + (1 << ofs_bits_s)}
  }, 1 << ofs_bits_l);

  const auto lz_memo = [&] {
    std::vector<std::array<encode::lz_data, 2>> ret(input.size());
    lz_helper lz_helper(input);
    for (size_t i = 0; i < input.size(); ++i) {
      encode::lz::find_all(i, ofs_tab, lz_min_len, ret[i], [&](size_t max_ofs) {
        return lz_helper.find(i, max_ofs, lz_min_len);
      });
      lz_helper.add_element(i);
    }
    return ret;
  }();

  auto sf_bitlens = [&]{
    // [Todo] Find a reasonable initialization.
    std::vector<size_t> ret(0x100 + lz_lens.size(), 18);
    const auto h = encode::huffman(utility::freq_u8(input), true);
    for (const auto w : h.words) ret[w] = h.code[w].bitlen;
    return ret;
  }();

  const auto update_costs = [&](const shannon_fano& sf, const size_t penalty = 1) {
    const auto& codewords = sf.codewords;
    size_t longest = 0;
    for (const auto w : sf.words) longest = std::max<size_t>(longest, codewords[w].bitlen);
    sf_bitlens.assign(sf_bitlens.size(), longest + penalty);
    for (const auto w : sf.words) sf_bitlens[w] = codewords[w].bitlen;
  };

  std::vector<uint8_t> best;

  while (true) {
    solver<tag> dp(input.size());
    for (size_t i = input.size(); i-- > 0; ) {
      dp.update_matrix(i, ofs_tab, lz_lens,
        [&](size_t li) { return sf_bitlens[lz_lens[li] + sf_offset]; },
        [&](size_t oi) { return lz_memo[i][oi]; },
        [&](size_t oi, size_t) -> tag { return {lz, oi}; }
      );
      dp.update(i, 1, sf_bitlens[input[i]], {uncomp, 0});
    }

    std::array<size_t, 0x100 + lz_lens.size()> code_counts = {};
    {
      size_t adr = 0;
      for (const auto& cmd : dp.optimal_path()) {
        switch (cmd.type.tag) {
        case uncomp: code_counts[input[adr]] += 1; break;
        case lz: code_counts[cmd.len + sf_offset] += 1; break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());
    }

    bool updated = false;
    for (size_t trial = 0; trial < 2; ++trial) {
      const bool add_header_cost = (trial == 0) ? false : true;
      auto sf = shannon_fano_encode(code_counts, add_header_cost);

      using namespace data_type;
      writer_b16_h ret;

      if (!dsp3) ret.write<bnh>({16, input.size()});
      ret.write<bnh>({16, sf.words.size()});
      ret.write<bnh>({16, dp.optimal_path().size()});

      size_t previous_word = 0x100 + lz_lens.size(); // (i.e. invalid word)
      for (const auto word : sf.words) {
        if (word <= previous_word || word >= previous_word + 0x14) {
          ret.write<bnh>({11, word});
        } else {
          size_t d = word - previous_word;
          if (d == 1) {
            ret.write<bnh>({2, 1});
          } else if (d <= 3) {
            ret.write<bnh>({3, 4 | (d - 2)});
          } else {
            ret.write<bnh>({6, 0x30 | (d - 4)});
          }
        }
        previous_word = word;
      }

      assert(sf.bits.size() == 4 || sf.bits.size() == 8);
      ret.write<b1>(sf.bits.size() == 4 ? 0 : 1);
      for (const auto b : sf.bits) ret.write<bnh>({3, b - 1});

      size_t adr = 0;
      for (const auto& cmd : dp.optimal_path()) {
        const auto [tag, oi] = cmd.type;
        switch (tag) {
        case uncomp: {
          ret.write<bnh>(sf.codewords[input[adr]]);
        } break;
        case lz: {
          ret.write<bnh>(sf.codewords[cmd.len + sf_offset]);
          const auto& o = ofs_tab[oi];
          ret.write<bnh>({o.bitlen, o.val + ((adr - cmd.lz_ofs()) - o.min)});
        } break;
        default: assert(0);
        }
        adr += cmd.len;
      }
      assert(adr == input.size());
      if (best.empty() || ret.size() < best.size()) {
        best = std::move(ret.out);
        update_costs(sf);
        updated = true;
        break;
      }
    }
    if (!updated) break;
  }

  return best;
}

std::vector<uint8_t> sd_gundam_gx_comp(std::span<const uint8_t> input) {
  return sd_gundam_gx_comp_core(input, true);
}

std::vector<uint8_t> slayers_comp(std::span<const uint8_t> input) {
  return sd_gundam_gx_comp_core(input, false);
}

} // namespace sfc_comp
