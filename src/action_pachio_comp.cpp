#include "lzss.hpp"

namespace sfc_comp {

std::vector<uint8_t> action_pachio_comp_core(
    std::span<const uint8_t> input, const size_t pad, const bool always_compress) {
  check_size(input.size(), 0, 0x8000);
  auto ret = lzss<writer_b8_l>(
    input, pad, [pad](std::span<uint8_t> in) {
      for (size_t i = 0; i < pad; ++i) in[i] = 0x20;
    },
    0x1000, 3, 0x12,
    4, true,
    [&] (size_t, size_t o, size_t l) {
      size_t d = (o - 0x12 - pad) & 0x0fff;
      return (d & 0x00ff) | (l - 3) << 8 | (d & 0x0f00) << 4;
    }
  );
  write16(ret, 0, input.size());
  if (input.size() == 0) {
    write16(ret, 2, 0xffff); // should be a (signed 16-bit) negative integer.
  } else {
    if (always_compress) {
      if (ret.size() - 4 >= 0x8000) {
        throw std::runtime_error("This algorithm cannot compress the given data.");
      }
    } else {
      if (ret.size() >= input.size()) {
        ret.resize(input.size() + 4);
        std::copy(input.begin(), input.end(), ret.begin() + 4);
      }
    }
    write16(ret, 2, ret.size() - 4);
  }
  return ret;
}

std::vector<uint8_t> super_dunk_star_comp(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x01, true);
}

std::vector<uint8_t> action_pachio_comp(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x04, true);
}

std::vector<uint8_t> keiba_eight_special_2_comp(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x06, true);
}

std::vector<uint8_t> dekitate_high_school_comp_1(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x0c, true);
}

std::vector<uint8_t> keirin_ou_comp(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x10, true);
}

std::vector<uint8_t> dekitate_high_school_comp_2(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x12, true);
}

std::vector<uint8_t> love_quest_comp(std::span<const uint8_t> input) {
  return action_pachio_comp_core(input, 0x12, false);
}

} // namespace sfc_comp
