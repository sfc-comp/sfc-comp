#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> kiki_kaikai_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = lzss<writer_b16_l>(
    input, 0, [](std::span<const uint8_t>) {},
    0x800, 3, 0x22,
    3, false,
    [&](size_t i, size_t o, size_t l) {
      return (l - 3) << 11 | ((i - o) - 1);
    }
  );
  ret[0] = 0x01;
  write16(ret, 1, input.size());
  if (input.size() > 0 && input.size() + 3 <= ret.size()) {
    ret.resize(input.size() + 3);
    std::ranges::copy(input, ret.begin() + 3);
    ret[0] = 0x00;
  }
  return ret;
}

namespace {

std::vector<uint8_t> wild_guns_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum tag { uncomp, rle, inc, dec };
  solver<tag> dp(input.size() / 2 * 2);
  auto c0_2 = dp.c<0, 2>(0x82); auto c2_2 = dp.c<2, 2>(0x80);

  size_t rlen16 = 0;
  for (size_t i = input.size() / 2 * 2; i > 0; ) {
    i -= 2;
    dp.update(i, 2, 0x80, c2_2, 1, uncomp);
    rlen16 = encode::run_length16_delta_r(input, i, rlen16);
    if (rlen16 >= 4) {
      const uint16_t delta = read16(input, i + 2) - read16(input, i);
      if (delta == 0x0001) dp.update(i, 4, 0x82, rlen16, c0_2, 3, inc);
      else if (delta == 0x0000) dp.update(i, 4, 0x82, rlen16, c0_2, 3, rle);
      else if (delta == 0xffff) dp.update(i, 4, 0x82, rlen16, c0_2, 3, dec);
    }
    c0_2.update(i); c2_2.update(i);
  }

  using namespace data_type;
  writer ret(3);
  if (input.size() & 1) ret.write<d8>(input.back());

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    switch (cmd.type) {
    case uncomp: ret.write<d8, d8n>(0x00 | (cmd.len / 2 - 1), {cmd.len, &input[adr]}); break;
    case rle: ret.write<d8, d16>(0x40 | (cmd.len / 2 - 2), read16(input, adr)); break;
    case inc: ret.write<d8, d16>(0x80 | (cmd.len / 2 - 2), read16(input, adr)); break;
    case dec: ret.write<d8, d16>(0xc0 | (cmd.len / 2 - 2), read16(input, adr)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret[0] = 0x02;
  write16(ret.out, 1, input.size());
  assert(adr == input.size() / 2 * 2);
  assert(dp.optimal_cost() + 3 + (input.size() & 1) == ret.size());
  return ret.out;
}

} // namespace

std::vector<uint8_t> wild_guns_comp(std::span<const uint8_t> input) {
  auto best = kiki_kaikai_comp(input);

  if (auto res = wild_guns_comp_2(input); best.empty() || res.size() < best.size()) {
    best = std::move(res);
  }

  return best;
}

} // namespace sfc_comp
