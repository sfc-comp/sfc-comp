#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

std::vector<uint8_t> kiki_kaikai_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum CompType { uncomp, lz };

  lz_helper lz_helper(input);
  sssp_solver<CompType> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), uncomp);
    auto res_lz = lz_helper.find_best(i, 0x800);
    dp.update_lz(i, 3, 0x22, res_lz, Constant<17>(), lz);
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b16 ret; ret.write<d8, d16>(0, 0);
  size_t adr = 0;
  for (const auto cmd : dp.commands()) {
    size_t d = adr - cmd.lz_ofs;
    switch (cmd.type) {
    case uncomp: ret.write<b1l, d8>(false, input[adr]); break;
    case lz: ret.write<b1l, d16>(true, (cmd.len - 3) << 11 | (d - 1)); break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.out[0] = 0x01;
  write16(ret.out, 1, input.size());
  assert((dp.total_cost() + 7) / 8 + 3 == ret.size() ||
         (dp.total_cost() + 7) / 8 + 4 == ret.size());
  assert(adr == input.size());

  if (input.size() > 0 && input.size() + 3 <= ret.size()) {
    ret.out.resize(input.size() + 3);
    std::copy(input.begin(), input.end(), ret.out.begin() + 3);
    ret[0] = 0x00;
  }
  return ret.out;
}

namespace {

std::vector<uint8_t> wild_guns_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);

  enum CompType { uncomp, rle, inc, dec };
  sssp_solver<CompType> dp(input.size() / 2 * 2);

  encode::rle_data rlen16 = {0, 0};
  for (size_t i = 0; i + 1 < input.size(); i += 2) {
    dp.update_k<2>(i, 2, 0x80, Linear<1, 1>(), uncomp);
    rlen16 = encode::run_length16_delta(input, i, rlen16);
    if (rlen16.v == 0x0001) dp.update_k<2>(i, 4, 0x82, rlen16.len, Constant<3>(), inc);
    else if (rlen16.v == 0x0000) dp.update_k<2>(i, 4, 0x82, rlen16.len, Constant<3>(), rle);
    else if (rlen16.v == 0xffff) dp.update_k<2>(i, 4, 0x82, rlen16.len, Constant<3>(), dec);
  }

  using namespace data_type;
  writer ret; ret.write<d8, d16>(0, 0);
  if (input.size() & 1) ret.write<d8>(input.back());

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
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
  assert(dp.total_cost() + 3 + (input.size() & 1) == ret.size());
  return ret.out;
}

} // namespace

std::vector<uint8_t> wild_guns_comp(std::span<const uint8_t> input) {
  auto best = kiki_kaikai_comp(input);

  if (auto res = wild_guns_comp_2(input); res.size() < best.size()) {
    best = std::move(res);
  }

  return best;
}

} // namespace sfc_comp
