#include <exception>
#include "MpiExit.h"
#include <mpi.h>

namespace optimet {
namespace {
bool &did_done_do_init() {
  static bool nrefs = false;
  return nrefs;
}

t_uint &mpi_reference() {
  static t_uint nrefs = 0;
  return nrefs;
}
} // anonymous namespace

void mpi_init(int argc, char **argv) {
  if(did_done_do_init())
    return;
  MPI_Init(&argc, &argv);
  did_done_do_init() = true;
}

bool mpi_initialized() { return did_done_do_init(); }

bool mpi_finalized() {
  int finalized;
  MPI_Finalized(&finalized);
  return finalized;
}

void mpi_finalize() {
  if(mpi_finalized() or not mpi_initialized())
    return;
  MPI_Finalize();
}

void increment_mpi_ref() { ++mpi_reference(); }
void decrement_mpi_ref() { --mpi_reference(); }
} /* optimet  */
