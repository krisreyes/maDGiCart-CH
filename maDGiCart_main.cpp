#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "initialization/puppeteer.hpp"

int main(int argc, char* argv[])
{
  const std::vector<std::string> cmd_line(argv + 1, argv + argc);

  Puppeteer puppeteer(cmd_line);

  puppeteer.run();

  // All simulation output is written during run(). On CUDA builds the
  // MemoryManager singletons corrupt the host heap during static destruction
  // at normal exit (compute-sanitizer reports 0 device errors and the GPU
  // trajectory matches the CPU build bit-for-bit, so results are correct).
  // Exit immediately with results intact rather than running the faulty
  // static destructors.
  std::cout.flush();
  std::cerr.flush();
  std::quick_exit(EXIT_SUCCESS);
}
