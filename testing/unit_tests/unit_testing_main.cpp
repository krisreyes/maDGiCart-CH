#include <cstdlib>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "logger/logger.hpp"
#include "testing/command_line.hpp"


int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // make command line arguments available to tests.
  std::vector<std::string> cmd_line(argv + 1, argv + argc);
  cmdline = cmd_line;

  const int rc = RUN_ALL_TESTS();

  // On CUDA builds the MemoryManager singletons corrupt the host heap during
  // static destruction at exit; all test output is already emitted, so exit
  // immediately with the gtest result code rather than running the destructors.
  std::cout.flush();
  std::cerr.flush();
  std::quick_exit(rc);
}

std::vector<std::string> cmdline;
