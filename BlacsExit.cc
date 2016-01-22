#include <exception>
#include <tuple>
#include "BlacsExit.h"
#include "Blacs.h"

namespace optimet {
namespace {
t_uint &blacs_reference() {
  static t_uint nrefs = 0;
  return nrefs;
}

std::tuple<t_uint, t_uint> blacs_pinfo() {
  static bool first_call = true;
  static std::tuple<t_uint, t_uint> result;
  if(first_call) {
    first_call = false;
    int rank, size;
    OPTIMET_FC_GLOBAL(blacs_pinfo, BLACS_PINFO)(&rank, &size);
    result = {static_cast<t_uint>(rank), static_cast<t_uint>(size)};
  }
  return result;
}
} // anonymous namespace
void blacks_exit(t_int status) {
  if(blacs_reference() != 0)
    throw std::runtime_error("At least one context has not been deleted");
  int stat = status;
  OPTIMET_FC_GLOBAL(blacs_exit, BLACS_EXIT)(&stat);
}

t_uint blacs_rank() { return std::get<0>(blacs_pinfo()); }
t_uint blacs_size() { return std::get<1>(blacs_pinfo()); }

void increment_blacs_ref() { ++blacs_reference(); }
void decrement_blacs_ref() { --blacs_reference(); }
} /* optimet  */
