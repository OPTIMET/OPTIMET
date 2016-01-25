#ifndef OPTIMET_MPI_TYPES_H
#define OPTIMET_MPI_TYPES_H

#include <mpi.h>
#include <complex>
#include "Types.h"

namespace optimet {
namespace mpi {
//! Type of an mpi tupe
typedef decltype(MPI_CHAR) MPIType;

//! MPI type associated with a c++ type
template <class T> struct Type;

#define OPTIMET_MACRO(TYPE) \
    template <> struct Type<TYPE> { static const MPIType value; };
  OPTIMET_MACRO(int8_t);
  OPTIMET_MACRO(int16_t);
  OPTIMET_MACRO(int32_t);
  OPTIMET_MACRO(int64_t);
  OPTIMET_MACRO(uint8_t);
  OPTIMET_MACRO(uint16_t);
  OPTIMET_MACRO(uint32_t);
  OPTIMET_MACRO(uint64_t);

  OPTIMET_MACRO(char);
  OPTIMET_MACRO(signed long);
  OPTIMET_MACRO(unsigned long);

  OPTIMET_MACRO(float);
  OPTIMET_MACRO(double);
  OPTIMET_MACRO(long double);
  OPTIMET_MACRO(std::complex<float>);
  OPTIMET_MACRO(std::complex<double>);
  OPTIMET_MACRO(std::complex<long double>);
#undef OPTIMET_MACRO

//! MPI type associated with a c++ type
template <class T> inline MPIType registered_type(T const&) { return Type<T>::value; }
} /* optime::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_TYPES */
