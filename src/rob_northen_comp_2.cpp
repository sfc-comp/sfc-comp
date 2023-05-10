#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> rob_northen_comp_2_core(
    std::span<const uint8_t> input,
    const size_t header_size, std::span<const vrange> ofs_tab) {
  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, uncompl, lz, lz2 };
  using tag = tag_ol<method>;

  static constexpr auto len_tab = to_vranges({
              // 0b10 (len == 2)
    {0x0003,  3, 0b110,          0b000},
    {0x0004,  3, 0b000,          0b010},
    {0x0006,  4, 0b0010,         0b0101},
    {0x0009, 11, 0b11100000001,  0b00011111111}
  }, 0x00ff);

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_k<4>(i, 12, 72, Linear<8, 9>(), {uncompl, 0, 0});
    const auto res_lz2 = lz_helper.find(i, 0x100, 2);
    dp.update_lz(i, 2, 2, res_lz2, Constant<11>(), {lz2, 0, 0});
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      1
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;
  writer_b8_h ret(header_size);
  ret.write<b1, b1>(false, false);

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    const size_t d = adr - cmd.lz_ofs;
    switch (cmd.type.tag) {
    case uncomp: {
      ret.write<b1, d8>(false, input[adr]);
    } break;
    case uncompl: {
      ret.write<bnh, d8n>({9, 0x170 + (cmd.len - 12) / 4}, {cmd.len, &input[adr]});
    } break;
    case lz2: {
      ret.write<b1, bnh, d8>(true, {2, 2}, d - 1);
    } break;
    case lz: {
      const auto& l = len_tab[cmd.type.li];
      const size_t ld = cmd.len - l.min;
      const auto lval = masked_add(l.val, ld, l.mask);
      ret.write<b1>(true);
      if (l.bitlen < 8) {
        ret.write<bnh>({l.bitlen, lval});
      } else {
        ret.write<bnh, d8>({l.bitlen - 8, l.val >> 8}, lval & 0xff);
      }
      const auto& o = ofs_tab[cmd.type.oi];
      const auto od = d - o.min;
      ret.write<bnh>({o.bitlen - 8, masked_add(o.val, od >> 8, o.mask)});
      ret.write<d8>(od & 0xff);
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  ret.write<bnh, d8, b1>({4, 0x0f}, 0, false);
  assert(header_size * 8 + 2 + dp.optimal_cost() + 4 + 8 + 1 == ret.bit_length());
  return ret.out;
}

static constexpr auto rnc2_offsets = to_vranges({
  {0x0001,  9, 0b0,      0b0},
  {0x0101, 11, 0b110,    0b000},
  {0x0201, 12, 0b1000,   0b0001},
  {0x0401, 13, 0b10101,  0b01010},
  {0x0801, 14, 0b101000, 0b010101},
}, 0x1000);

} // namespace

std::vector<uint8_t> rob_northen_comp_2(std::span<const uint8_t> input) {
  auto ret = rob_northen_comp_2_core(input, 0x12, rnc2_offsets);
  write32b(ret, 0, 0x524e4302);
  write32b(ret, 4, input.size());
  write32b(ret, 8, ret.size() - 18);
  write16b(ret, 12, utility::crc16(input, 0, input.size()));
  write16b(ret, 14, utility::crc16(ret, 18, ret.size() - 18));
  ret[0x10] = 0; // leeway
  ret[0x11] = 0; // chunks
  return ret;
}

std::vector<uint8_t> spirou_comp(std::span<const uint8_t> input) {
  auto ret = rob_northen_comp_2_core(input, 0x02, rnc2_offsets);
  write16b(ret, 0, 0x524e);
  return ret;
}

std::vector<uint8_t> smurfs_comp(std::span<const uint8_t> input) {
  auto ret = rob_northen_comp_2_core(input, 0x02, std::span(rnc2_offsets.begin(), 3));
  write16b(ret, 0, 0x524e);
  return ret;
}

} // namespace sfc_comp
