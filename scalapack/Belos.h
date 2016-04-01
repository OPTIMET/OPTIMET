#ifndef OPTIMET_SCALAPACK_BELOS_H_
#define OPTIMET_SCALAPACK_BELOS_H_

#include "Types.h"
#ifdef OPTIMET_BELOS

#include "mpi/Communicator.h"
#include "scalapack/Matrix.h"

#include <Kokkos_Serial.hpp>
#include <Kokkos_View.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Tpetra_Vector.hpp>

namespace optimet {
namespace scalapack {

//! Simplifies access to a unmanaged kokkos view type
template <class SCALAR, class DEVICE = Kokkos::Serial>
using TeuchosArrayView = Teuchos::ArrayView<SCALAR>;
//! Simplifies access to a tpetra vector type
template <class SCALAR> using TpetraVector = Tpetra::Vector<SCALAR>;

//! Obtains a 1-dimensional kokkos view over a matrix data
template <class SCALAR> TeuchosArrayView<SCALAR> view(Matrix<SCALAR> &A);
//! Obtains a 1-dimensional kokkos view over a constant matrix data
template <class SCALAR> TeuchosArrayView<SCALAR const> view(Matrix<SCALAR> const &A);

//! \brief Obtains a map for a scalapack::Matrix
//! \note Belos doesn't know that the matrix is block cyclic. Since we will be defining our own
//! matrix multiply, that's not really necessary.
//! \param[in] A: matrix
//! \param[in] comm: an mpi communicator equivalent to the scalapack context. All processes and only
//!     the processes for which the scalapack context is valid should be in the communicator.
template <class SCALAR>
Teuchos::RCP<const Tpetra::Map<>>
matrix_map(Matrix<SCALAR> const &A, mpi::Communicator const &comm);

//! \brief Obtains a vector for a scalapack matrix
//! \param[in] A: A matrix viewed as a row vector.
//! \param[in] comm: an mpi communicator equivalent to the scalapack context. All processes and only
//!     the processes for which the scalapack context is valid should be in the communicator.
template <class SCALAR>
TpetraVector<SCALAR> tpetra_vector(Matrix<SCALAR> const &A, mpi::Communicator const &comm);
}
}

#include "scalapack/belos.hpp"
#endif
#endif
