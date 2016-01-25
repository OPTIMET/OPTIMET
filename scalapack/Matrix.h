#ifndef OPTIMET_SCALAPACK_MATRIX_H_
#define OPTIMET_SCALAPACK_MATRIX_H_

#include <array>
#include "Types.h"
#include "scalapack/Context.h"

namespace optimet {
namespace scalapack {
class Matrix {
public:
  //! Underlying scalar type
  typedef t_real Scalar;
  //! Underlying eigen matrix type
  typedef optimet::Matrix<Scalar> EigenMatrix;
  //! Rows and Colums of the local blocks
  struct Sizes {
    t_uint rows, cols;
  };
  //! Indices of the process starting the distribution
  struct Index {
    t_uint row, col;
  };

  //! Constructs from any eigen matrix
  Matrix(Context const &context, Sizes size, Sizes blocks, Index index = {0, 0})
      : context_(context), blacs_{{1, context.is_valid() ? *context : -1,
                                   static_cast<int>(size.rows), static_cast<int>(size.cols),
                                   static_cast<int>(blocks.rows), static_cast<int>(blocks.cols),
                                   static_cast<int>(index.row), static_cast<int>(index.col)}},
        matrix_(EigenMatrix::Zero(rows(context, size, blocks, index),
                                  cols(context, size, blocks, index))) {
    blacs_[8] = EigenMatrix::IsRowMajor ? eigen().rows() : eigen().cols();
  };

  //! Gets the underlying eigen matrix
  EigenMatrix &eigen() { return matrix_; }
  //! Gets the underlying eigen matrix
  EigenMatrix const &eigen() const { return matrix_; }

  //! The underlying blacs context
  Context const &context() { return context_; }

  //! Gets the underlying blacs construct
  std::array<int, 9> const &blacs() const { return blacs_; }

  //! Transfer matrix to another context
  // Matrix transfer(Context other) const;

protected:
  //! Associated blacs context
  Context context_;
  //! Blacs construct
  std::array<int, 9> blacs_;
  //! The underlying eigen matrix
  EigenMatrix matrix_;

  static EigenMatrix::Index rows(Context const &context, Sizes size, Sizes blocks, Index index);
  static EigenMatrix::Index cols(Context const &context, Sizes size, Sizes blocks, Index index);
};
}
}
#endif
