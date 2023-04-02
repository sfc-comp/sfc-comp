#include "io.hpp"
#include "utility.hpp"

namespace sfc_comp {

namespace io {

std::vector<uint8_t> load(const std::string& path) {
  std::ifstream is(path, std::ios::binary);
  if (!is.is_open()) {
    throw std::runtime_error(format("Cannot open \"%s\".", path.c_str()));
  }
  is.seekg(0, is.end);
  size_t length = is.tellg();
  is.seekg(0, is.beg);
  std::vector<uint8_t> ret(length);
  is.read(reinterpret_cast<char*>(ret.data()), length);
  is.close();
  return ret;
}

void save(const std::string& path, std::span<const uint8_t> data) {
  std::ofstream os(path, std::ios::binary);
  if (!os.is_open()) {
    throw std::runtime_error(format("Cannot open \"%s\".", path.c_str()));
  }
  os.write(reinterpret_cast<const char*>(data.data()), data.size());
  os.close();
}

} // namespace input

} // namespace sfc_comp
