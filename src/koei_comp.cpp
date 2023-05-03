#include "algorithm.hpp"
#include "encode.hpp"
#include "utility.hpp"
#include "writer.hpp"

namespace sfc_comp {

namespace {

template <typename Writer8, typename WriterRaw, typename Writer16>
void koei_comp_core(std::span<const uint8_t> input,
    Writer8& writer8, WriterRaw& writer_raw, Writer16& writer16) {

  check_size(input.size(), 0, 0x800000);

  enum method { uncomp, lz };
  using tag = tag_ol<method>;

  static constexpr auto ofs_tab = std::to_array<vrange>({
    vrange(0x0001, 0x0004,  6, 0b0000'00),
    vrange(0x0005, 0x0008,  7, 0b00010'00),
    vrange(0x0009, 0x0020,  8, 0b00011'000),
    vrange(0x0021, 0x0080,  9, 0b0011'00000),
    vrange(0x0081, 0x0100, 10, 0b011'0000000),
    vrange(0x0101, 0x0200, 11, 0b100'00000000),
    vrange(0x0201, 0x0400, 12, 0b101'000000000),
    vrange(0x0401, 0x0800, 13, 0b110'0000000000),
    vrange(0x0801, 0x1000, 14, 0b111'00000000000)
  });

  static constexpr auto len_tab = std::to_array<vrange>({
    vrange(0x0002, 0x0002,  1, 0b1),
    vrange(0x0003, 0x0004,  3, 0b01'0),
    vrange(0x0005, 0x0008,  5, 0b001'00),
    vrange(0x0009, 0x0010,  7, 0b0001'000),
    vrange(0x0011, 0x0020,  9, 0b00001'0000),
    vrange(0x0021, 0x0040, 11, 0b000001'00000),
    vrange(0x0041, 0x0080, 13, 0b0000001'000000),
    vrange(0x0081, 0x00ff, 14, 0b0000000'0000000)
  });

  lz_helper lz_helper(input);
  sssp_solver<tag> dp(input.size());

  for (size_t i = 0; i < input.size(); ++i) {
    dp.update(i, 1, 1, Constant<9>(), {uncomp, 0, 0});
    dp.update_lz_matrix(i, ofs_tab, len_tab,
      [&](size_t oi) { return lz_helper.find_best(i, ofs_tab[oi].max); },
      [&](size_t oi, size_t li) -> tag { return {lz, oi, li}; },
      1
    );
    lz_helper.add_element(i);
  }

  using namespace data_type;

  size_t adr = 0;
  for (const auto& cmd : dp.commands()) {
    switch (cmd.type.tag) {
    case uncomp: {
      writer8.template write<b1>(true);
      writer_raw.template write<d8>(input[adr]);
    } break;
    case lz: {
      const auto& l = len_tab[cmd.type.li];
      const auto& o = ofs_tab[cmd.type.oi];
      writer8.template write<b1>(false);
      writer16.template write<bnh>({l.bitlen, l.val + (cmd.len - l.min)});
      writer16.template write<bnh>({o.bitlen, o.val + ((adr - cmd.lz_ofs) - o.min)});
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
