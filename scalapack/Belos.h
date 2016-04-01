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

template <class SCALAR, class DEVICE = Kokkos::Serial>
using KokkosView = Kokkos::View<SCALAR *, DEVICE, Kokkos::LayoutLeft, Kokkos::MemoryUnmanaged>;

//! \brief Obtains a kokkos view over a matrix
//! \param[in] A: A matrix which must be a *row vector*
template <class SCALAR> KokkosView<SCALAR> view(Matrix<SCALAR> &A);
//! Obtains a kokkos view over a constant matrix
//! \param[in] A: A matrix which must be a *row vector*
template <class SCALAR> KokkosView<SCALAR const> view(Matrix<SCALAR> const &A);

//! \brief Obtains a map for a scalapack::Matrix
//! \note Belos doesn't know that the matrix is block cyclic. Since we will be defining our own
//! matrix multiply, that's not really necessary.
template <class SCALAR>
Teuchos::RCP<const Tpetra::Map<>> matrix_map(Matrix<SCALAR> const &, mpi::Communicator const &);
}
}

#include "scalapack/belos.hpp"
#endif
#endif
