#define CATCH_CONFIG_RUNNER

#include <catch.hpp>
#include <random>
#include <memory>
#include <mpi.h>
#include "MpiExit.h"

std::unique_ptr<std::mt19937_64> mersenne(new std::mt19937_64(0));

int main( int argc, char* const argv[] )
{
  Catch::Session session; // There must be exactly once instance

  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
    return returnCode;
  mersenne.reset(new std::mt19937_64(session.configData().rngSeed));

  optimet::mpi_init(argc, argv);

  auto const result = session.run();
  optimet::mpi_finalize();
  return result;
}
