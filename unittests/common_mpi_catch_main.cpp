#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_CONSOLE_WIDTH 100

#include "Types.h"
#include "mpi/Session.h"
#include <catch.hpp>
#include <memory>
#include <random>

std::unique_ptr<std::mt19937_64> mersenne(new std::mt19937_64(0));

int main(int argc, const char **argv) {
  Catch::Session session; // There must be exactly once instance

  int returnCode = session.applyCommandLine(argc, argv);
  if(returnCode != 0) // Indicates a command line error
    return returnCode;
  mersenne.reset(new std::mt19937_64(session.configData().rngSeed));

  optimet::mpi::init(argc, argv);

  auto const result = session.run();
  optimet::mpi::finalize();
  return result;
}
