#ifndef OPTIMET_SCALAPACK_COLLECTIVES_H
#define OPTIMET_SCALAPACK_COLLECTIVES_H
#include "Types.h"
#ifdef OPTIMET_MPI

#include <mpi.h>
#include <type_traits>
#include <vector>

namespace optimet {
namespace scalapack {
namespace details {
//! True if a type scalapack knows about
template <class T>
struct is_fundamental
    : std::integral_constant<bool, std::is_same<T, int>::value or std::is_same<T, float>::value or
                                       std::is_same<T, double>::value or
                                       std::is_same<T, std::complex<float>>::value or
                                       std::is_same<T, std::complex<double>>::value> {};
}

class Context;
#define OPTIMET_MACRO(TYPE)                                                                        \
  /** Broadcasts from given process to all others **/                                              \
  TYPE broadcast(TYPE const &, Context const &context, t_uint row, t_uint col);                    \
  Matrix<TYPE> broadcast(Matrix<TYPE> const &, Context const &context, t_uint row, t_uint col);    \
  Vector<TYPE> broadcast(Vector<TYPE> const &, Context const &context, t_uint row, t_uint col);
OPTIMET_MACRO(int);
OPTIMET_MACRO(float);
OPTIMET_MACRO(double);
OPTIMET_MACRO(std::complex<double>);
OPTIMET_MACRO(std::complex<float>);
#undef OPTIMET_MACRO
} /* optime::mpi */
} /* optimet */
#endif
#endif /* ifndef OPTIMET_MPI_COMMUNICATOR */
