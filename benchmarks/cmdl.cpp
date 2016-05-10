#include "Types.h"
#include "cmdl.h"

#ifdef OPTIMET_BELOS
#include <Teuchos_CommandLineProcessor.hpp>

namespace optimet {
namespace {
struct Cmdl {
  std::string solver = "scalapack";
  t_real tolerance = 1e-8;
  t_int maximum_restarts = 20;
  t_int itermax = 1000;
  t_int block_size = 1;
  t_int adaptive_block_size = 1;
  t_int verbosity = 0;
  t_int num_blocks = 50;
  t_int num_recycled_blocks = 5;
};
}
Teuchos::RCP<Teuchos::ParameterList> parse_cmdl(int argc, char *argv[]) {
  Teuchos::CommandLineProcessor clp(false, false); // don't throw exceptions
  Cmdl cmdl;

  clp.setOption("solver", &cmdl.solver, "Name of the solver to use");
  clp.setOption("tolerance", &cmdl.tolerance, "Convergence criteria");
  clp.setOption("maximum_restarts", &cmdl.maximum_restarts,
                "Number of times to restart the solver, if supported");
  clp.setOption("itermax", &cmdl.itermax, "Maximum number of iterations");
  clp.setOption("bock_size", &cmdl.block_size);
  clp.setOption("adaptive_bock_size", &cmdl.adaptive_block_size);
  clp.setOption("verbosity", &cmdl.verbosity);
  clp.setOption("num_blocks", &cmdl.num_blocks);
  clp.setOption("num_recycled_blocks", &cmdl.num_recycled_blocks);
  auto result = Teuchos::rcp(new Teuchos::ParameterList);
  if(clp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    throw std::runtime_error("Couldn't parse command-line");

  result->set("Solver", cmdl.solver);
  result->set("Convergence Tolerance", cmdl.tolerance);
  result->set("Maximum Restarts", cmdl.maximum_restarts);
  result->set("Maximum Iterations", cmdl.itermax);
  result->set("Num Blocks", cmdl.num_blocks);
  result->set("Num Recycled Blocks", cmdl.num_recycled_blocks);
  result->set("Block Size", cmdl.block_size);
  result->set("Adaptive Block Size", cmdl.adaptive_block_size);
  result->set("Verbosity", cmdl.verbosity);
  return result;
}
}
#endif
