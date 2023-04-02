#include <chrono>

#include "sfc_comp.hpp"
#include "version.h"

#define STRINGIFY(n) #n
#define TOSTRING(n) STRINGIFY(n)

int main(int argc, char** argv) {

#if defined(COMP_FUNC) && defined(COMP_EXT)
  using namespace sfc_comp;
  using namespace std::chrono;

  printf("%s (%s)\n\n", VERSION_STRING, TOSTRING(COMP_FUNC));
  for (int i = 1; i < argc; ++i) {
    const std::string src_path = argv[i];
    const auto dest_path = src_path + TOSTRING(COMP_EXT);

    try {
      const auto beg = high_resolution_clock::now();

      const auto input = io::load(src_path);
      const auto output = COMP_FUNC(input);
      io::save(dest_path, output);

      const auto end = high_resolution_clock::now();

      printf("   %s\n-> %s\n", src_path.c_str(), dest_path.c_str());
      printf("     elapsed time: %.3f sec\n", duration_cast<nanoseconds>(end - beg).count() / 1e9);
      printf("uncompressed size: %4zXh Byte(s)\n", input.size());
      printf("  compressed size: %4zXh Byte(s)\n", output.size());
      if (input.size() > 0) {
        printf("compression ratio: %5.2f%%\n", output.size() * 100. / input.size());
      }
    } catch (const std::exception& e) {
      printf("[Error] Failed to compress: %s\n", src_path.c_str());
      printf("// %s\n", e.what());
    }
    puts("");
  }
#endif

  printf("Press Enter to exit.\n");
  (void) getchar();

  return 0;
}
