#pragma once

#include <cstddef>
#include <cstdint>

#include <fstream>
#include <vector>

#include <span>

namespace sfc_comp {

namespace data_type {

struct none {

};

struct b1l {
  b1l(bool b) : b(b) {}
  bool b;
};

struct b1h {
  b1h(bool b) : b(b) {}
  bool b;
};

struct b4h {
  b4h(uint8_t x) : x(x) {}
  uint8_t x;
};

struct b8ln_h {
  b8ln_h(size_t n, size_t x) : n(n), x(x) {}
  size_t n, x;
};

struct b8ln_l {
  b8ln_l(size_t n, size_t x) : n(n), x(x) {}
  size_t n, x;
};

struct b8hn_h {
  b8hn_h(size_t n, size_t x) : n(n), x(x) {}
  size_t n, x;
};

struct b8hn_l {
  b8hn_l(size_t n, size_t x) : n(n), x(x) {}
  size_t n, x;
};

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
  d24(uint32_t x) : x(x) {}
  uint32_t x;
};

struct d24b {
  d24b(uint32_t x) : x(x) {}
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
  writer() {}

  template <typename Head, typename... Args>
  void write(const Head& h, Args... args) {
    write(h);
    write(args...);
  }

  void write() {}

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

struct writer_b4 : public writer {
  writer_b4() : lo(false) {}

  template <typename Head, typename... Args>
  void write(const Head& h, Args... args) {
    write(h);
    write(args...);
  }

  void write() {}

  void write(const data_type::b4h& d) {
    if (!lo) {
      out.push_back(d.x << 4);
    } else {
      out.back() |= d.x;
    }
    lo = !lo;
  }

  void write(const data_type::d8& d) {
    write<data_type::b4h>(d.x >> 4);
    write<data_type::b4h>(d.x & 0x0f);
  }

  void write(const data_type::d16b& d) {
    write<data_type::d8>(d.x >> 8);
    write<data_type::d8>(d.x >> 0);
  }

  void write(const data_type::d8n& d) {
    for (size_t i = 0; i < d.n; ++i) write<data_type::d8>(d.ptr[i]);
  }

  bool lo;
};

struct writer_b : public writer {
  using writer::write;

  writer_b() : bit(0), bits_pos(0) {}

  template <typename Head, typename... Args>
  void write(const Head& h, Args... args) {
    write(h);
    write(args...);
  }

  void write(const data_type::b1l& d) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0);
      bit = 8;
    }
    --bit;
    if (d.b) out[bits_pos] |= 1 << (7 - bit);
  }

  void write(const data_type::b1h& d) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0);
      bit = 8;
    }
    --bit;
    if (d.b) out[bits_pos] |= 1 << bit;
  }

  void write(const data_type::none&) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0);
      bit = 8;
    }
  }

  void write(const data_type::b8ln_l& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1l>((d.x >> i) & 1);
    }
  }

  void write(const data_type::b8ln_h& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1l>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  void write(const data_type::b8hn_l& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1h>((d.x >> i) & 1);
    }
  }

  void write(const data_type::b8hn_h& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1h>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  size_t bit;
  size_t bits_pos;
};

struct writer_b16 : public writer {
  using writer::write;

  writer_b16() : bit(0), bits_pos(0) {}

  template <typename Head, typename... Args>
  void write(const Head& h, Args... args) {
    write(h);
    write(args...);
  }

  void write(const data_type::b1l& d) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0); out.push_back(0);
      bit = 16;
    }
    --bit;
    if (d.b) out[bits_pos + ((15 - bit) >> 3)] |= 1 << (7 - (bit & 7));
  }

  void write(const data_type::b1h& d) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0); out.push_back(0);
      bit = 16;
    }
    --bit;
    if (d.b) out[bits_pos + (bit >> 3)] |= 1 << (bit & 7);
  }

  void write(const data_type::none&) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0); out.push_back(0);
      bit = 16;
    }
  }

  void write(const data_type::b8ln_l& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1l>((d.x >> i) & 1);
    }
  }

  void write(const data_type::b8ln_h& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1l>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  void write(const data_type::b8hn_l& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1h>((d.x >> i) & 1);
    }
  }

  void write(const data_type::b8hn_h& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1h>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  size_t bit;
  size_t bits_pos;
};

struct writer_b16_hasty : public writer {
  using writer::write;

  writer_b16_hasty() : bit(0), bits_pos(0) {}

  template <typename Head, typename... Args>
  void write(const Head& h, Args... args) {
    write(h);
    write(args...);
  }

  void trim() {
    if (bit == 16 && bits_pos + 2 == out.size()) {
      out.pop_back(); out.pop_back();
      bit = 0;
    }
  }

  void write(const data_type::b1l& d) {
    --bit;
    if (d.b) out[bits_pos + ((15 - bit) >> 3)] |= 1 << (7 - (bit & 7));
    write(data_type::none());
  }

  void write(const data_type::b1h& d) {
    --bit;
    if (d.b) out[bits_pos + (bit >> 3)] |= 1 << (bit & 7);
    write(data_type::none());
  }

  void write(const data_type::none&) {
    if (bit == 0) {
      bits_pos = out.size();
      out.push_back(0); out.push_back(0);
      bit = 16;
    }
  }

  void write(const data_type::b8ln_l& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1l>((d.x >> i) & 1);
    }
  }

  void write(const data_type::b8ln_h& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1l>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  void write(const data_type::b8hn_l& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1h>((d.x >> i) & 1);
    }
  }

  void write(const data_type::b8hn_h& d) {
    for (size_t i = 0; i < d.n; ++i) {
      write<data_type::b1h>((d.x >> (d.n - 1 - i)) & 1);
    }
  }

  size_t bit;
  size_t bits_pos;
};

} // namespace sfc_comp
