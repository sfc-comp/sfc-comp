#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

template <typename Writer8, typename WriterRaw, typename Writer16>
requires std::derived_from<Writer8, writer> &&
         std::derived_from<WriterRaw, writer> &&
         std::derived_from<Writer16, writer>
void koei_comp_core(std::span<const uint8_t> input,
    Writer8& writer8, WriterRaw& writer_raw, Writer16& writer16) {

  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = to_vranges({
    {0x0001,  6, 0b0000'00},
    {0x0005,  7, 0b00010'00},
    {0x0009,  8, 0b00011'000},
    {0x0021,  9, 0b0011'00000},
    {0x0081, 10, 0b011'0000000},
    {0x0101, 11, 0b100'00000000},
    {0x0201, 12, 0b101'000000000},
    {0x0401, 13, 0b110'0000000000},
    {0x0801, 14, 0b111'00000000000}
  }, 0x1000);

  static constexpr auto len_tab = to_vranges({
    {0x0002,  1, 0b1},
    {0x0003,  3, 0b01'0},
    {0x0005,  5, 0b001'00},
    {0x0009,  7, 0b0001'000},
    {0x0011,  9, 0b00001'0000},
    {0x0021, 11, 0b000001'00000},
    {0x0041, 13, 0b0000001'000000},
    {0x0081, 14, 0b0000000'0000000}
  }, 0x00ff);

  lz_helper lz_helper(input, true);
  solver<tag> dp(input.size());
  auto c0 = dp.template c<0>(len_tab.back().max);

  for (size_t i = input.size(); i--; ) {
    lz_helper.reset(i);
    dp.update(i, 1, 9, {uncomp, 0, 0});
    dp.update_matrix(i, ofs_tab, len_tab, c0, 1,
      [&](size_t oi) { return lz_helper.find(i, ofs_tab[oi].max, len_tab.front().min); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; }
    );
    c0.update(i);
  }

  using namespace data_type;

  size_t adr = 0;
  for (const auto& cmd : dp.optimal_path()) {
    const auto [tag, oi, li] = cmd.type;
    switch (tag) {
    case uncomp: {
      writer8.template write<b1>(true);
      writer_raw.template write<d8>(input[adr]);
    } break;
    case lz: {
      const auto& l = len_tab[li];
      const auto& o = ofs_tab[oi];
      writer8.template write<b1>(false);
      writer16.template write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      writer16.template write<bnh>({o.bitlen, o.val + ((adr - cmd.lz_ofs()) - o.min)});
    } break;
    default: assert(0);
    }
    adr += cmd.len;
  }
  writer8.template write<b1>(false);
  writer16.template write<bnh>({14, 0b00000001111111});
  assert(adr == input.size());
}

struct writer_koei : public writer_b8_h {
  using writer_b8_h::write_;

  writer_koei(size_t s = 2) : writer_b8_h(s), bit16(0), curr16(-1), next16(0) {}

  template <typename Head, typename... Args>
  void write(const Head& h, const Args&... args) {
    write_(h);
    if constexpr (sizeof... (args) >= 1) {
      write(args...);
    }
  }

  void truncate() {
    if (next16 + 2 == size()) {
      out.pop_back(); out.pop_back();
    }
  }

 protected:
  void write_(const data_type::bnh& d) {
    for (size_t b = d.n; b-- > 0; ) {
      if (bit16 == 0) {
        curr16 = next16; next16 = size();
        out.push_back(0); out.push_back(0);
        bit16 = 16;
      }
      --bit16;
      if ((d.v >> b) & 1) {
        (*this)[curr16 + (bit16 >> 3)] |= 1 << (bit16 & 7);
      }
    }
  }

 public:
  size_t bit16;
  size_t curr16, next16;
};

} // namespace

std::vector<uint8_t> koei_comp(std::span<const uint8_t> input) {
  writer_koei ret;
  koei_comp_core(input, ret, ret, ret);
  ret.truncate();
  return ret.out;
}

std::vector<uint8_t> battle_cross_comp(std::span<const uint8_t> input) {
  writer_b8_h ret(4); writer raw; writer_b16_h flags;
  koei_comp_core(input, ret, raw, flags);
  write16(ret.out, 0, ret.size());
  ret.extend(raw);
  write16(ret.out, 2, ret.size() - 2);
  ret.extend(flags);
  return ret.out;
}

} // namespace sfc_comp
