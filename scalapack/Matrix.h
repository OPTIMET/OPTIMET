#ifndef OPTIMET_SCALAPACK_MATRIX_H_
#define OPTIMET_SCALAPACK_MATRIX_H_

#include "Types.h"
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include <array>

namespace optimet {
namespace scalapack {

#ifdef OPTIMET_MPI
//! Wrapper around scalapack distributed matrix and eigen
template <class SCALAR = t_real> class Matrix {
public:
  //! Underlying scalar type
  typedef SCALAR Scalar;
  //! Underlying eigen matrix type
  typedef optimet::Matrix<Scalar> EigenMatrix;
  //! Constructs from any eigen matrix
  Matrix(Context const &context, Sizes sizes, Sizes blocks, Index index = {0, 0})
      : context_(context),

        matrix_(EigenMatrix::Zero(rows(context, sizes, blocks, index),
                                  cols(context, sizes, blocks, index))),

        blacs_{{1, context.is_valid() ? *context : -1, static_cast<int>(sizes.rows),
                static_cast<int>(sizes.cols), static_cast<int>(blocks.rows),
                static_cast<int>(blocks.cols), static_cast<int>(index.row),
                static_cast<int>(index.col), static_cast<int>(local_leading())}} {}
  //! Copy constructor
  Matrix(Matrix const &other)
      : context_(other.context_), matrix_(other.matrix_), blacs_(other.blacs_) {}
  //! Move Copy constructor
  Matrix(Matrix &&other)
      : context_(std::move(other.context_)), matrix_(std::move(other.matrix_)),
        blacs_(std::move(other.blacs_)) {}

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

  //! \brief Assigns matrix
  //! \details Copies from other to current. The two matrices must have the same size.
  //! They may have different contexts.
  void operator=(Matrix<SCALAR> const &other);
  //! Move assignment
  void operator=(Matrix<SCALAR> &&other);

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

  //! Process row index over which the first row of the matrix is distributed
  t_uint first_row() const { return static_cast<t_uint>(blacs_[6]); }
  //! Process column index over which the first column of the matrix is distributed
  t_uint first_col() const { return static_cast<t_uint>(blacs_[7]); }
};

//! Multiplies two matrices
template <class SCALAR>
void pdgemm(SCALAR alpha, Matrix<SCALAR> const &a, Matrix<SCALAR> const &b, SCALAR beta,
            Matrix<SCALAR> &c, char opa = 'N', char opb = 'N');

#endif
}
}

#ifdef OPTIMET_MPI
#include "scalapack/Matrix.hpp"
#endif
#endif
