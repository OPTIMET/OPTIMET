#ifndef OPTIMET_SCALAPACK_BELOS_HPP_
#define OPTIMET_SCALAPACK_BELOS_HPP_

#include "Types.h"
#ifdef OPTIMET_BELOS

#include "scalapack/Belos.h"
#include <Kokkos_Serial.hpp>
#include <Kokkos_View.hpp>

namespace optimet {
namespace scalapack {
template <class SCALAR> TeuchosArrayView<SCALAR> view(Matrix<SCALAR> &A) {
  return TeuchosArrayView<SCALAR>(A.local().data(), A.local().size());
}

template <class SCALAR> TeuchosArrayView<SCALAR const> view(Matrix<SCALAR> const &A) {
  return TeuchosArrayView<SCALAR const>(A.local().data(), A.local().size());
}

template <class SCALAR>
Teuchos::RCP<const Tpetra::Map<>>
matrix_map(Matrix<SCALAR> const &A, mpi::Communicator const &comm) {
  typedef int ordinal;

  // Construct teuchos mpi wrapper, making sure it owns it.
  MPI_Comm raw_comm;
  MPI_Comm_dup(*comm, &raw_comm);
  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
      Teuchos::opaqueWrapper(raw_comm, MPI_Comm_free);
  Teuchos::RCP<const Teuchos::Comm<ordinal>> tcomm =
      Teuchos::rcp(new Teuchos::MpiComm<ordinal>(opaque_comm));

  // Now creates map
  Teuchos::RCP<const Tpetra::Map<>> result =
      Teuchos::rcp(new Tpetra::Map<>(A.size(), A.local().size(), 0, tcomm));

  // And populate it
  return result;
}

template <class SCALAR>
TpetraVector<SCALAR> tpetra_vector(Matrix<SCALAR> const &A, mpi::Communicator const &comm) {
  return TpetraVector<SCALAR>(matrix_map(A, comm), view<SCALAR>(A));
}
}
}
#endif
#endif
