#include <exception>
#include <tuple>
#include "scalapack/InitExit.h"
#include "scalapack/Blacs.h"

namespace optimet {
namespace scalapack {
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
    OPTIMET_FC_GLOBAL_(blacs_pinfo, BLACS_PINFO)(&rank, &size);
    result = std::tuple<t_uint, t_uint>{static_cast<t_uint>(rank), static_cast<t_uint>(size)};
  }
  return result;
}

bool &has_did_done_exit() {
  static bool has = false;
  return has;
}

} // anonymous namespace

bool finalized() { return has_did_done_exit(); }
void finalize(t_int status) {
  if(finalized())
    return;
  has_did_done_exit() = true;
  int stat = status;
  OPTIMET_FC_GLOBAL_(blacs_exit, BLACS_EXIT)(&stat);
}

t_uint global_rank() { return std::get<0>(blacs_pinfo()); }
t_uint global_size() { return std::get<1>(blacs_pinfo()); }

void increment_ref() { ++blacs_reference(); }
void decrement_ref() { --blacs_reference(); }
} /* scalapack  */
} /* optimet  */
