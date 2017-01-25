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

#ifndef OPTIMET_SCALAPACK_MATRIX_H_
#define OPTIMET_SCALAPACK_MATRIX_H_

#include "Types.h"

#ifdef OPTIMET_SCALAPACK
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include <array>

namespace optimet {
namespace scalapack {

//! \brief Traits associated with Matrix
//! \details Allows specialization for Eigen::Map
template <class SCALAR> struct MatrixTraits {
  //! Underlying scalar type
  typedef SCALAR Scalar;
  //! Underlying type of the Eigen matrix
  typedef optimet::Matrix<Scalar> EigenMatrix;
};
//! Specialization when memory is held outside Eigen
template <class SCALAR> struct MatrixTraits<SCALAR *> {
  //! Underlying scalar type
  typedef SCALAR Scalar;
  //! Underlying type of the Eigen matrix
  typedef Eigen::Map<typename MatrixTraits<Scalar>::EigenMatrix> EigenMatrix;
};
//! Specialization when memory is held outside Eigen
template <class SCALAR> struct MatrixTraits<SCALAR const *> {
  //! Underlying scalar type
  typedef SCALAR Scalar;
  //! Underlying type of the Eigen matrix
  typedef Eigen::Map<const typename MatrixTraits<Scalar>::EigenMatrix> EigenMatrix;
};

//! Wrapper around scalapack distributed matrix and eigen
template <class SCALAR = t_real> class Matrix {
public:
  //! Traits associated with this object
  typedef MatrixTraits<SCALAR> traits;
  //! Underlying scalar type
  typedef typename traits::Scalar Scalar;
  //! Underlying eigen matrix type
  typedef typename traits::EigenMatrix EigenMatrix;
  //! Matrix where memory is held by Eigen
  typedef Matrix<Scalar> ConcreteMatrix;
  //! Matrix where memory is held outside Eigen
  typedef Matrix<Scalar *> MapMatrix;
  //! Copy constructor
  Matrix(Matrix const &other)
      : context_(other.context_), matrix_(other.matrix_), blacs_(other.blacs_) {}
  //! Move Copy constructor
  Matrix(Matrix &&other)
      : context_(std::move(other.context_)), matrix_(std::move(other.matrix_)),
        blacs_(std::move(other.blacs_)) {}

  Matrix(EigenMatrix const &matrix, Context const &context, Sizes sizes, Sizes blocks,
         Index index = {0, 0})
      : context_(context), matrix_(matrix),
        blacs_{{1, context.is_valid() ? *context : -1, static_cast<int>(sizes.rows),
                static_cast<int>(sizes.cols), static_cast<int>(blocks.rows),
                static_cast<int>(blocks.cols), static_cast<int>(index.row),
                static_cast<int>(index.col), static_cast<int>(local_leading())}} {
    assert(blocks.rows > 0);
    assert(blocks.cols > 0);
#ifndef NDEBUG
    auto const nrows = local_rows(context, sizes, blocks, index);
    auto const ncols = local_cols(context, sizes, blocks, index);
    assert(nrows == matrix_.rows());
    assert(ncols == matrix_.cols());
#endif
  }
  //! Constructs from any eigen matrix
  Matrix(Context const &context, Sizes sizes, Sizes blocks, Index index = {0, 0})
      : Matrix(EigenMatrix::Zero(local_rows(context, sizes, blocks, index),
                                 local_cols(context, sizes, blocks, index)),
               context, sizes, blocks, index) {}

  //! Gets the underlying local eigen matrix
  EigenMatrix &local() { return matrix_; }
  //! Gets the underlying local eigen matrix
  EigenMatrix const &local() const { return matrix_; }

  //! The underlying blacs context
  Context const &context() const { return context_; }

  //! Gets the underlying blacs construct
  std::array<int, 9> const &blacs() const { return blacs_; }

  //! Transfer matrix to another context
  void transfer_to(Context const &un, ConcreteMatrix &other) const;
  //! Transfer matrix to another context
  void transfer_to(ConcreteMatrix &other) const {
    return transfer_to(context().size() > other.size() ? context() : other.context(), other);
  }
  //! Transfer matrix to another context
  void transfer_to(Context const &un, MapMatrix &other) const;
  //! Transfer matrix to another context
  void transfer_to(MapMatrix &other) const {
    return transfer_to(context().size() > other.size() ? context() : other.context(), other);
  }
  //! Transfer matrix to another context
  ConcreteMatrix transfer_to(Context const &un, Context const &other,
                             Sizes const &blocs = {std::numeric_limits<t_uint>::max(),
                                                   std::numeric_limits<t_uint>::max()},
                             Index const &index = {std::numeric_limits<t_uint>::max(),
                                                   std::numeric_limits<t_uint>::max()}) const;
  ConcreteMatrix transfer_to(Context const &other,
                             Sizes const &blocks = {std::numeric_limits<t_uint>::max(),
                                                    std::numeric_limits<t_uint>::max()},
                             Index const &index = {std::numeric_limits<t_uint>::max(),
                                                   std::numeric_limits<t_uint>::max()}) const {
    return transfer_to(context().size() > other.size() ? context() : other, other, blocks, index);
  }
  //! \details Copies from other to current. The two matrices must have the same size.
  //! They may have different contexts.
  void operator=(Matrix const &other);
  //! Move assignment
  void operator=(Matrix &&other);

  //! Index of the process with th upper left element
  Index index() const { return {static_cast<t_uint>(blacs()[6]), static_cast<t_uint>(blacs()[7])}; }
  //! Size of the cyclic blocks
  Sizes blocks() const {
    return {static_cast<t_uint>(blacs()[4]), static_cast<t_uint>(blacs()[5])};
  }
  //! Size of the distributed matrix
  Sizes sizes() const { return {static_cast<t_uint>(blacs()[2]), static_cast<t_uint>(blacs()[3])}; }
  //! Number of rows of the distributed matrix
  t_uint rows() const { return static_cast<t_uint>(blacs()[2]); }
  //! Number of rows of the distributed matrix
  t_uint cols() const { return static_cast<t_uint>(blacs()[3]); }
  //! Number of rows of the distributed matrix
  t_uint size() const { return rows() * cols(); }

  //! Local leading dimension
  t_uint local_leading() const {
    auto const result = EigenMatrix::IsRowMajor ? local().cols() : local().rows();
    return std::max(result, static_cast<decltype(result)>(1));
  }

  //! Global to local indices
  std::tuple<t_uint, t_uint, t_uint, t_uint> local_indices(t_uint i, t_uint j) const {
    return local_indices(std::make_tuple(i, j));
  }
  //! Global to local indices
  std::tuple<t_uint, t_uint, t_uint, t_uint>
  local_indices(std::tuple<t_uint, t_uint> const &i) const;

  //! Local to global indices
  //! \param[in] i: local row index
  //! \param[in] j: local column index
  //! \param[in] proc_i: process row index
  //! \param[in] proc_j: process column index
  std::tuple<t_uint, t_uint>
  global_indices(t_uint i, t_uint j, t_uint proc_i, t_uint proc_j) const {
    return global_indices(std::make_tuple(i, j, proc_i, proc_j));
  }
  std::tuple<t_uint, t_uint> global_indices(t_uint i, t_uint j) const {
    return global_indices(std::make_tuple(i, j, context().row(), context().col()));
  }
  // Local to global indices
  std::tuple<t_uint, t_uint>
  global_indices(std::tuple<t_uint, t_uint, t_uint, t_uint> const &i) const;

  //! Sizes of local matrix
  Sizes local_sizes() const {
    return {static_cast<t_uint>(local().rows()), static_cast<t_uint>(local().cols())};
  }
  //! Computes local size for given global size
  Sizes local_sizes(Sizes const &gsizes) const {
    return {static_cast<t_uint>(local_rows(context(), gsizes, blocks(), index())),
            static_cast<t_uint>(local_cols(context(), gsizes, blocks(), index()))};
  }

  //! Number of local rows for given scalapack parameters
  static typename EigenMatrix::Index
  local_rows(Context const &context, Sizes size, Sizes blocks, Index index);
  //! Number of local columns for given scalapack parameters
  static typename EigenMatrix::Index
  local_cols(Context const &context, Sizes size, Sizes blocks, Index index);

protected:
  //! Associated blacs context
  Context context_;
  //! The underlying eigen matrix
  EigenMatrix matrix_;
  //! Blacs construct
  std::array<int, 9> blacs_;
  //! Process row index over which the first row of the matrix is distributed
  t_uint first_row() const { return static_cast<t_uint>(blacs_[6]); }
  //! Process column index over which the first column of the matrix is distributed
  t_uint first_col() const { return static_cast<t_uint>(blacs_[7]); }
};

//! Creates a scalapack map of input matrix
template <class SCALAR>
Matrix<SCALAR *> map_matrix(optimet::Matrix<SCALAR> &matrix, Context const &context,
                            Sizes const &sizes, Sizes const &blocks);
//! Creates a scalapack map of input matrix
template <class SCALAR>
Matrix<SCALAR const *> map_cmatrix(optimet::Matrix<SCALAR> const &matrix, Context const &context,
                                   Sizes const &sizes, Sizes const &blocks) {
  return map_matrix(matrix, context, sizes, blocks);
}
//! Creates a scalapack map of input matrix
template <class SCALAR>
Matrix<SCALAR const *> map_matrix(optimet::Matrix<SCALAR> const &matrix, Context const &context,
                                  Sizes const &sizes, Sizes const &blocks);
//! Creates a scalapack map of input vector
template <class SCALAR>
Matrix<SCALAR *> map_matrix(optimet::Vector<SCALAR> &matrix, Context const &context,
                            Sizes const &sizes, Sizes const &blocks);
//! Creates a scalapack map of input vector
template <class SCALAR>
Matrix<SCALAR const *> map_matrix(optimet::Vector<SCALAR> const &matrix, Context const &context,
                                  Sizes const &sizes, Sizes const &blocks);
//! Creates a scalapack map of input matrix
template <class SCALAR>
Matrix<SCALAR const *> map_cmatrix(optimet::Vector<SCALAR> const &matrix, Context const &context,
                                   Sizes const &sizes, Sizes const &blocks) {
  return map_matrix(matrix, context, sizes, blocks);
}

//! Multiplies two matrices
template <class SCALAR>
void pdgemm(typename MatrixTraits<SCALAR>::Scalar alpha, Matrix<SCALAR> const &a,
            Matrix<SCALAR> const &b, typename MatrixTraits<SCALAR>::Scalar beta, Matrix<SCALAR> &c,
            char opa = 'N', char opb = 'N');
//! Multiplies two matrices
template <class SCALAR>
void pdgemm(typename MatrixTraits<SCALAR>::Scalar alpha, Matrix<SCALAR> const &a,
            Matrix<SCALAR const *> const &b, typename MatrixTraits<SCALAR>::Scalar beta,
            Matrix<SCALAR *> &c, char opa = 'N', char opb = 'N');
//! Multiplies two matrices
template <class SCALAR>
void pdgemm(typename MatrixTraits<SCALAR>::Scalar alpha, Matrix<SCALAR const *> const &a,
            Matrix<SCALAR const *> const &b, typename MatrixTraits<SCALAR>::Scalar beta,
            Matrix<SCALAR *> &c, char opa = 'N', char opb = 'N');
}
}

#include "scalapack/Matrix.hpp"
#endif
#endif
