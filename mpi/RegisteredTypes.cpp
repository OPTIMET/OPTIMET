// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "mpi/RegisteredTypes.h"

namespace optimet {
namespace mpi {
MPIType const Type<std::int8_t>::value = MPI_INT8_T;
MPIType const Type<std::int16_t>::value = MPI_INT16_T;
MPIType const Type<std::int32_t>::value = MPI_INT32_T;
MPIType const Type<std::int64_t>::value = MPI_INT64_T;
MPIType const Type<std::uint8_t>::value = MPI_UINT8_T;
MPIType const Type<std::uint16_t>::value = MPI_UINT16_T;
MPIType const Type<std::uint32_t>::value = MPI_UINT32_T;
MPIType const Type<std::uint64_t>::value = MPI_UINT64_T;

#ifndef OPTIMET_CHAR_ARCH
MPIType const Type<char>::value = MPI_CHAR;
#endif
#ifndef OPTIMET_LONG_ARCH
MPIType const Type<signed long>::value = MPI_LONG;
#endif
#ifndef OPTIMET_ULONG_ARCH
MPIType const Type<unsigned long>::value = MPI_UNSIGNED_LONG;
#endif

MPIType const Type<float>::value = MPI_FLOAT;
MPIType const Type<double>::value = MPI_DOUBLE;
MPIType const Type<long double>::value = MPI_LONG_DOUBLE;
MPIType const Type<std::complex<float>>::value = MPI_C_FLOAT_COMPLEX;
MPIType const Type<std::complex<double>>::value = MPI_C_DOUBLE_COMPLEX;
MPIType const Type<std::complex<long double>>::value = MPI_C_LONG_DOUBLE_COMPLEX;

static_assert(is_registered_type<int>::value, "Checking int is registered");
static_assert(is_registered_type<std::complex<double>>::value,
              "Checking complex double is registered");
static_assert(not is_registered_type<std::complex<int>>::value,
              "Checking complex int is NOT registered");
} /* mpi */
} /* optimet */
