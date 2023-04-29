#include "lzss.hpp"

namespace sfc_comp {

namespace {

std::vector<uint8_t> dq12_comp_core(std::span<const uint8_t> input, const size_t header_size, const size_t max_dist) {
  auto ret = lzss<writer_b8_l>(
    input, 0x12, [](std::span<const uint8_t>) {},
    max_dist, 3, 0x12,
    header_size, true,
    [&](size_t, size_t o, size_t l) {
      size_t d = (o - 0x24) & 0xfff;
      return (d & 0x00ff) | (l - 3) << 8 | (d & 0x0f00) << 4;
    }
  );
  return ret;
}

} // namespace

std::vector<uint8_t> lemmings_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  auto ret = dq12_comp_core(input, 2, 0x1000);
  write16(ret, 0, ret.size() - 2);
  return ret;
}

std::vector<uint8_t> maka_maka_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x8000);
  auto ret = dq12_comp_core(input, 2, 0x0fff);
  write16(ret, 0, ret.size() - 2);
  return ret;
}

std::vector<uint8_t> gionbana_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x8000);
  auto ret = dq12_comp_core(input, 3, 0x1000);
  if (ret.size() - 3 > 0x8000) {
    throw std::runtime_error("This algorithm cannot compress the given data.");
  }
  write24(ret, 0, ret.size() - 3); // unknown
  return ret;
}

std::vector<uint8_t> final_stretch_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0xffff);
  auto ret = dq12_comp_core(input, 4, 0x1000);
  write16(ret, 0, input.size());
  if (input.size() > 0) {
    write16(ret, 2, ret.size() - 4);
  } else {
    write16(ret, 2, 0x0001);
  }
  return ret;
}

std::vector<uint8_t> dq12_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = dq12_comp_core(input, 2, 0x1000);
  write16b(ret, 0, input.size());
  return ret;
}

std::vector<uint8_t> dq5_comp_2(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = dq12_comp_core(input, 2, 0x1000);
  write16(ret, 0, input.size());
  return ret;
}

std::vector<uint8_t> odekake_lester_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0xffff);
  auto ret = dq12_comp_core(input, 2, 0x0fff); // cf. $CE:81C4, $CE:8210
  write16(ret, 0, input.size());
  return ret;
}

std::vector<uint8_t> ihatovo_monogatari_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 1, 0x10000);
  auto ret = dq12_comp_core(input, 4, 0x1000);
  write32(ret, 0, input.size()); // unknown
  return ret;
}

std::vector<uint8_t> super_jinsei_game_comp(std::span<const uint8_t> input) {
  check_size(input.size(), 0, 0x10000);
  return dq12_comp_core(input, 0, 0x1000);
}

} // namespace sfc_comp
