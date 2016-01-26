#ifndef OPTIMET_SCALAPACK_MATRIX_H_
#define OPTIMET_SCALAPACK_MATRIX_H_

#include <array>
#include "Types.h"
#include "scalapack/Context.h"

namespace optimet {
namespace scalapack {

//! Rows and Colums of the local blocks
struct Sizes {
  t_uint rows, cols;
};
//! Indices of the process starting the distribution
struct Index {
  t_uint row, col;
};

#ifdef OPTIMET_MPI
//! Wrapper around scalapack distributed matrix and eigen
template <class SCALAR = t_real> class Matrix {
public:
  //! Underlying scalar type
  typedef SCALAR Scalar;
  //! Underlying eigen matrix type
  typedef optimet::Matrix<Scalar> EigenMatrix;
  //! Constructs from any eigen matrix
  Matrix(Context const &context, Sizes size, Sizes blocks, Index index = {0, 0})
      : context_(context),

        matrix_(EigenMatrix::Zero(rows(context, size, blocks, index),
                                  cols(context, size, blocks, index))),

        blacs_{{1, context.is_valid() ? *context : -1, static_cast<int>(size.rows),
                static_cast<int>(size.cols), static_cast<int>(blocks.rows),
                static_cast<int>(blocks.cols), static_cast<int>(index.row),
                static_cast<int>(index.col), static_cast<int>(local_leading())}} {}

  //! Gets the underlying local eigen matrix
  EigenMatrix &local() { return matrix_; }
  //! Gets the underlying local eigen matrix
  EigenMatrix const &local() const { return matrix_; }

  //! The underlying blacs context
  Context const &context() const { return context_; }

  //! Gets the underlying blacs construct
  std::array<int, 9> const &blacs() const { return blacs_; }

  //! Transfer matrix to another context
  void transfer_to(Context const &un, Matrix &other) const;
  //! Transfer matrix to another context
  void transfer_to(Matrix &other) const {
    return transfer_to(context().size() > other.size() ? context() : other.context(), other);
  }
  //! Transfer matrix to another context
  Matrix transfer_to(Context const &un, Context const &other,
                     Sizes const &blocs = {std::numeric_limits<t_uint>::max(),
                                           std::numeric_limits<t_uint>::max()},
                     Index const &index = {std::numeric_limits<t_uint>::max(),
                                           std::numeric_limits<t_uint>::max()}) const;
  Matrix transfer_to(Context const &other,
                     Sizes const &blocks = {std::numeric_limits<t_uint>::max(),
                                            std::numeric_limits<t_uint>::max()},
                     Index const &index = {std::numeric_limits<t_uint>::max(),
                                           std::numeric_limits<t_uint>::max()}) const {
    return transfer_to(context().size() > other.size() ? context() : other, other, blocks, index);
  }

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
  t_uint local_leading() const { return EigenMatrix::IsRowMajor ? local().cols() : local().rows(); }

protected:
  //! Associated blacs context
  Context context_;
  //! The underlying eigen matrix
  EigenMatrix matrix_;
  //! Blacs construct
  std::array<int, 9> blacs_;

  static typename EigenMatrix::Index
  rows(Context const &context, Sizes size, Sizes blocks, Index index);
  static typename EigenMatrix::Index
  cols(Context const &context, Sizes size, Sizes blocks, Index index);
};
#endif
}
}

#ifdef OPTIMET_MPI
#include "scalapack/Matrix.hpp"
#endif
#endif
