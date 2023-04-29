#pragma once

#include <cstddef>
#include <cstdint>

#include <fstream>
#include <vector>

#include <span>

namespace sfc_comp {

namespace data_type {

struct d8 {
  d8(uint8_t x) : x(x) {}
  uint8_t x;
};

struct d16 {
  d16(uint16_t x) : x(x) {}
  uint16_t x;
};

struct d16b {
  d16b(uint16_t x) : x(x) {}
  uint16_t x;
};

struct d24 {
  d24(uint32_t x) : x(x & 0xffffff) {}
  uint32_t x;
};

struct d24b {
  d24b(uint32_t x) : x(x & 0xffffff) {}
  uint32_t x;
};

struct d32 {
  d32(uint32_t x) : x(x) {}
  uint32_t x;
};

struct d32b {
  d32b(uint32_t x) : x(x) {}
  uint32_t x;
};

struct d8n {
  d8n(size_t n, const uint8_t* ptr) : n(n), ptr(ptr) {}
  size_t n;
  const uint8_t* ptr;
};

struct d8n_l {
  d8n_l(size_t n, const uint8_t* ptr) : n(n), ptr(ptr) {}
  size_t n;
  const uint8_t* ptr;
};

struct d8n_h {
  d8n_h(size_t n, const uint8_t* ptr) : n(n), ptr(ptr) {}
  size_t n;
  const uint8_t* ptr;
};

struct d8n_lb {
  d8n_lb(size_t n, const uint8_t* ptr) : n(n), ptr(ptr) {}
  size_t n;
  const uint8_t* ptr;
};

struct d8n_hb {
  d8n_hb(size_t n, const uint8_t* ptr) : n(n), ptr(ptr) {}
  size_t n;
  const uint8_t* ptr;
};

struct d8nk {
  d8nk(size_t n, const uint8_t* ptr, size_t delta) :
    n(n), ptr(ptr), delta(delta) {}
  size_t n;
  const uint8_t* ptr;
  size_t delta;
};

} // namespace data_type

struct writer {
  writer(size_t s = 0) : out(s) {}

  template <typename Head, typename... Args>
  void write(const Head& h, const Args&... args) {
    write(h);
    if constexpr (sizeof... (args) >= 1) {
      write(args...);
    }
  }

  void write(const data_type::d8& d) {
    out.push_back(d.x);
  }

  void write(const data_type::d16& d) {
    out.push_back(d.x);
    out.push_back(d.x >> 8);
  }

  void write(const data_type::d16b& d) {
    out.push_back(d.x >> 8);
    out.push_back(d.x);
  }

  void write(const data_type::d24& d) {
    out.push_back(d.x);
    out.push_back(d.x >> 8);
    out.push_back(d.x >> 16);
  }

  void write(const data_type::d24b& d) {
    out.push_back(d.x >> 16);
    out.push_back(d.x >> 8);
    out.push_back(d.x);
  }

  void write(const data_type::d32& d) {
    out.push_back(d.x);
    out.push_back(d.x >> 8);
    out.push_back(d.x >> 16);
    out.push_back(d.x >> 24);
  }

  void write(const data_type::d32b& d) {
    out.push_back(d.x >> 24);
    out.push_back(d.x >> 16);
    out.push_back(d.x >> 8);
    out.push_back(d.x);
  }

  void write(const data_type::d8n& d) {
    for (size_t i = 0; i < d.n; ++i) out.push_back(d.ptr[i]);
  }

  void write(const data_type::d8nk& d) {
    for (size_t i = 0; i < d.n; i += d.delta) out.push_back(d.ptr[i]);
  }

  void write(const data_type::d8n_lb& d) {
    for (size_t i = 0; i < (d.n & ~1); i += 2) {
      out.push_back((d.ptr[i] << 4) | (d.ptr[i + 1] & 0x0f));
    }
    if (d.n & 1) out.push_back(d.ptr[d.n - 1] << 4);
  }

  void write(const data_type::d8n_hb& d) {
    for (size_t i = 0; i < (d.n & ~1); i += 2) {
      out.push_back((d.ptr[i] & 0xf0) | (d.ptr[i + 1] >> 4));
    }
    if (d.n & 1) out.push_back(d.ptr[d.n - 1] & 0xf0);
  }

  uint8_t& operator [] (size_t i) {
    return out[i];
  }

  size_t size() const {
    return out.size();
  }

  std::vector<uint8_t> out;
};

namespace data_type {

struct none {

};

struct b1 {
  b1(bool b) : b(b) {}
  bool b;
};

struct b4 {
  b4(uint8_t x) : x(x & 0x0f) {}
  uint8_t x;
};

struct bnh {
  bnh(size_t n, size_t x) : n(n), x(x) {}
  size_t n, x;
};

struct bnl {
  bnl(size_t n, size_t x) : n(n), x(x) {}
  size_t n, x;
};

} // namespace data_type

template <bool LowNibbleFirst>
struct writer_b4 : public writer {
  writer_b4(size_t s = 0) : writer(s), nibble(0), nibble_pos(-1) {}

  template <typename Head, typename... Args>
  void write(const Head& h, const Args&... args) {
    write(h);
    if constexpr (sizeof... (args) >= 1) {
      write(args...);
    }
  }

  void write(const data_type::none&) {
    if (nibble == 0) {
      nibble_pos = out.size();
      out.push_back(0);
      nibble = 2;
    }
  }

  void write(const data_type::b4& d) {
    write<data_type::none>({});
    --nibble;
    if constexpr (LowNibbleFirst) {
      out[nibble_pos] |= (d.x & 0x0f) << (4 * (1 - nibble));
    } else {
      out[nibble_pos] |= (d.x & 0x0f) << (4 * nibble);
    }
  }

  void write(const data_type::d8& d) {
    write<data_type::b4>(d.x >> 4);
    write<data_type::b4>(d.x & 0x0f);
  }

  void write(const data_type::d16b& d) {
    write<data_type::d8>(d.x >> 8);
    write<data_type::d8>(d.x >> 0);
  }

  void write(const data_type::d8n& d) {
    for (size_t i = 0; i < d.n; ++i) write<data_type::d8>(d.ptr[i]);
  }

  size_t nibble_size() const {
    return size() * 2 - nibble;
  }

  size_t nibble;
  size_t nibble_pos;
};

using writer_b4_l = writer_b4<true>;
using writer_b4_h = writer_b4<false>;

template <size_t BlockBytes, bool LSBFirst, bool PreRead>
struct bitstream_writer : public writer {
  static constexpr size_t block_bitsize = BlockBytes * 8;
  using writer::write;

  bitstream_writer(size_t s = 0) : writer(s), bit(0), bits_pos(-1) {}

  template <typename Head, typename... Args>
  void write(const Head& h, const Args&... args) {
    write(h);
    if constexpr (sizeof... (args) >= 1) {
      write(args...);
    }
  }

  void write(const data_type::b1& d) {
    if constexpr (!PreRead) write(data_type::none());
    --bit;
    if constexpr (LSBFirst) {
      const size_t b = block_bitsize - 1 - bit;
      if (d.b) out[bits_pos + b / 8] |= 1 << (b % 8);
    } else {
      if (d.b) out[bits_pos + bit / 8] |= 1 << (bit % 8);
    }
    if constexpr (PreRead) write(data_type::none());
  }

  void write(const data_type::none&) {
    if (bit == 0) {
      bits_pos = out.size();
      for (size_t i = 0; i < BlockBytes; ++i) out.push_back(0);
      bit = block_bitsize;
    }
  }

  void write(const data_type::bnl& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1>((d.x >> i) & 1);
    }
  }

  void write(const data_type::bnh& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  void trim() {
    if (bit == block_bitsize && bits_pos + BlockBytes == out.size()) {
      for (size_t i = 0; i < BlockBytes; ++i) out.pop_back();
      bit = 0;
    }
  }

  size_t bit_length() const {
    return size() * 8 - bit;
  }

  size_t bit;
  size_t bits_pos;
};

template <bool LSBFirst> using writer_b = bitstream_writer<1, LSBFirst, false>;
using writer_b_l = writer_b<true>;
using writer_b_h = writer_b<false>;

template <bool LSBFirst> using writer_b8 = bitstream_writer<1, LSBFirst, false>;
using writer_b8_l = writer_b8<true>;
using writer_b8_h = writer_b8<false>;

template <bool LSBFirst> using writer_b16 = bitstream_writer<2, LSBFirst, false>;
using writer_b16_l = writer_b16<true>;
using writer_b16_h = writer_b16<false>;

template <bool LSBFirst> using writer_b16_hasty = bitstream_writer<2, LSBFirst, true>;
using writer_b16_hasty_l = writer_b16_hasty<true>;
using writer_b16_hasty_h = writer_b16_hasty<false>;

} // namespace sfc_comp
