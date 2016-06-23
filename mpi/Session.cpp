#include "mpi/Session.h"
#include <exception>
#include <mpi.h>
#ifdef OPTIMET_BELOS
#include <Tpetra_Core.hpp>
#endif

namespace optimet {
namespace mpi {
namespace {
bool &did_done_do_init() {
  static bool nrefs = false;
  return nrefs;
}

t_uint &global_reference() {
  static t_uint nrefs = 0;
  return nrefs;
}
} // anonymous namespace

void init(int argc, const char **argv) {
  if(did_done_do_init())
    return;
#ifdef OPTIMET_BELOS
  Tpetra::initialize(&argc, const_cast<char ***>(&argv));
#else
  MPI_Init(&argc, const_cast<char ***>(&argv));
#endif
  did_done_do_init() = true;
}

bool initialized() { return did_done_do_init(); }

bool finalized() {
  int finalized;
  MPI_Finalized(&finalized);
  return finalized;
}

void finalize() {
  if(finalized() or not initialized())
    return;
#ifdef OPTIMET_BELOS
  Kokkos::finalize();
# endif
  MPI_Finalize();
}

void increment_ref() { ++global_reference(); }
void decrement_ref() { --global_reference(); }
} /* optimet::mpi  */
} /* optimet  */
