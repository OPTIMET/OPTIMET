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

static_assert(not std::is_same<char, std::int8_t>::value, "");
#define OPTIMET_MACRO(TYPE) \
    template <> struct Type<TYPE> { static const MPIType value; };
  OPTIMET_MACRO(std::int8_t);
  OPTIMET_MACRO(std::int16_t);
  OPTIMET_MACRO(std::int32_t);
  OPTIMET_MACRO(std::int64_t);
  OPTIMET_MACRO(std::uint8_t);
  OPTIMET_MACRO(std::uint16_t);
  OPTIMET_MACRO(std::uint32_t);
  OPTIMET_MACRO(std::uint64_t);

#ifndef OPTIMET_CHAR_ARCH
  OPTIMET_MACRO(char);
#endif
#ifndef OPTIMET_LONG_ARCH
  OPTIMET_MACRO(signed long);
#endif
#ifndef OPTIMET_ULONG_ARCH
  OPTIMET_MACRO(unsigned long);
#endif

  OPTIMET_MACRO(float);
  OPTIMET_MACRO(double);
  OPTIMET_MACRO(long double);
  OPTIMET_MACRO(std::complex<float>);
  OPTIMET_MACRO(std::complex<double>);
  OPTIMET_MACRO(std::complex<long double>);
#undef OPTIMET_MACRO

//! MPI type associated with a c++ type
template <class T> inline MPIType registered_type(T const&) { return Type<T>::value; }

namespace details {
//! Defines c++17 metafunction
template <class...> using void_t = void;
}
//! True if the type is registered
template <class T, class = details::void_t<>> class is_registered_type : public std::false_type {};
template <class T>
class is_registered_type<T, details::void_t<decltype(Type<T>::value)>> : public std::true_type {};

static_assert(is_registered_type<int>::value, "Checking int is registered");
static_assert(is_registered_type<std::complex<double>>::value,
              "Checking complex double is registered");
static_assert(not is_registered_type<std::complex<int>>::value,
              "Checking complex int is NOT registered");
} /* optime::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_TYPES */
