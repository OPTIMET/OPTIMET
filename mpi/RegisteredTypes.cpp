#include "mpi/RegisteredTypes.h"

namespace optimet {
namespace mpi {
MPIType const Type<int8_t>::value   = MPI_INT8_T;
MPIType const Type<int16_t>::value  = MPI_INT16_T;
MPIType const Type<int32_t>::value  = MPI_INT32_T;
MPIType const Type<int64_t>::value  = MPI_INT64_T;
MPIType const Type<uint8_t>::value   = MPI_UINT8_T;
MPIType const Type<uint16_t>::value  = MPI_UINT16_T;
MPIType const Type<uint32_t>::value  = MPI_UINT32_T;
MPIType const Type<uint64_t>::value  = MPI_UINT64_T;

MPIType const Type<char>::value           = MPI_CHAR;
MPIType const Type<signed long>::value    = MPI_LONG;
MPIType const Type<unsigned long>::value  = MPI_UNSIGNED_LONG;

MPIType const Type<float>::value        = MPI_FLOAT;
MPIType const Type<double>::value       = MPI_DOUBLE;
MPIType const Type<long double>::value  = MPI_LONG_DOUBLE;
MPIType const Type<std::complex<float>>::value       = MPI_C_FLOAT_COMPLEX;
MPIType const Type<std::complex<double>>::value      = MPI_C_DOUBLE_COMPLEX;
MPIType const Type<std::complex<long double>>::value = MPI_C_LONG_DOUBLE_COMPLEX;
} /* mpi */
} /* optimet */
