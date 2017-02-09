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

#ifndef OPTIMET_SCALAPACK_BELOS_HPP_
#define OPTIMET_SCALAPACK_BELOS_HPP_

#ifndef OPTIMET_SCALAPACK_BELOS_H_
#error This file should not be included explicitly. Use Belos.h instead.
#endif

#include "Types.h"
#ifdef OPTIMET_BELOS

#include "scalapack/Belos.h"
#include "scalapack/Matrix.h"
#include <Kokkos_Serial.hpp>
#include <Kokkos_View.hpp>

namespace optimet {
namespace scalapack {
template <class SCALAR> TeuchosArrayView<typename Matrix<SCALAR>::Scalar> view(Matrix<SCALAR> &A) {
  typedef typename Matrix<SCALAR>::Scalar Scalar;
  return TeuchosArrayView<Scalar>(A.local().data(), A.local().size());
}

template <class SCALAR>
TeuchosArrayView<typename Matrix<SCALAR>::Scalar const> view(Matrix<SCALAR> const &A) {
  typedef typename Matrix<SCALAR>::Scalar Scalar;
  return TeuchosArrayView<Scalar const>(A.local().data(), A.local().size());
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
Teuchos::RCP<TpetraVector<typename Matrix<SCALAR>::Scalar>>
tpetra_vector(Matrix<SCALAR> const &A, mpi::Communicator const &comm) {
  typedef typename Matrix<SCALAR>::Scalar Scalar;
  return Teuchos::rcp(
      new TpetraVector<Scalar>(matrix_map(A, comm), view<SCALAR>(A), A.local().size(), 1));
}

template <class SCALAR>
void matrix_vector_operator(Matrix<SCALAR> const &A,
                            const TpetraVector<typename Matrix<SCALAR>::Scalar> &X,
                            TpetraVector<typename Matrix<SCALAR>::Scalar> &Y, Belos::ETrans trans) {
  if(not A.context().is_valid())
    return;
  if(X.getLocalLength() != Y.getLocalLength())
    throw std::runtime_error("Local lengths of X and Y are different");
  char character = 'N';
  if(trans == Belos::TRANS)
    character = 'T';
  else if(trans == Belos::CONJTRANS)
    character = 'C';
  auto output_map = as_matrix(Y, A);
  auto const input_map = as_matrix(X, A);
  pdgemm(1e0, A, input_map, 0e0, output_map, character, 'N');
}

template <class SCALAR>
Matrix<SCALAR *> as_matrix(TpetraVector<SCALAR> &X, Context const &context, Sizes const &blocks,
                           Index const &index) {
  auto const data = X.getLocalLength() ? X.getDataNonConst(0).getRawPtr() : nullptr;
  Sizes const global_sizes{static_cast<t_uint>(X.getGlobalLength()),
                           static_cast<t_uint>(X.getNumVectors())};
  Sizes const local_sizes{
      static_cast<t_uint>(Matrix<SCALAR>::local_rows(context, global_sizes, blocks, index)),
      static_cast<t_uint>(Matrix<SCALAR>::local_cols(context, global_sizes, blocks, index))};
  typedef typename Matrix<SCALAR>::EigenMatrix Map;
  return {Map::Map(data, local_sizes.rows, local_sizes.cols), context, global_sizes, blocks, index};
}

template <class SCALAR>
Matrix<SCALAR const *> as_matrix(TpetraVector<SCALAR> const &X, Context const &context,
                                 Sizes const &blocks, Index const &index) {
  auto const data = X.getLocalLength() ? X.getData(0).getRawPtr() : nullptr;
  Sizes const global_sizes{static_cast<t_uint>(X.getGlobalLength()),
                           static_cast<t_uint>(X.getNumVectors())};
  Sizes const local_sizes{
      static_cast<t_uint>(Matrix<SCALAR>::local_rows(context, global_sizes, blocks, index)),
      static_cast<t_uint>(Matrix<SCALAR>::local_cols(context, global_sizes, blocks, index))};
  typedef typename Matrix<SCALAR>::EigenMatrix Map;
  return {Map::Map(data, local_sizes.rows, local_sizes.cols), context, global_sizes, blocks, index};
}
}
}
#endif
#endif
