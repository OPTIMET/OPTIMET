#include <exception>
#include "BlacsExit.h"
#include "Blacs.h"

namespace optimet {
namespace {
t_uint &blacs_reference() {
  static t_uint nrefs = 0;
  return nrefs;
}
} // anonymous namespace
void blacks_exit(t_int status) {
  if(blacs_reference() != 0)
    throw std::runtime_error("At least one context has not been deleted");
  int stat = status;
  OPTIMET_FC_GLOBAL(blacs_exit, BLACS_EXIT)(&stat);
}

void increment_blacs_ref() { ++blacs_reference(); }
void decrement_blacs_ref() { --blacs_reference(); }
} /* optimet  */
