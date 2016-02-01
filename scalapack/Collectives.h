#ifndef OPTIMET_SCALAPACK_COLLECTIVES_H
#define OPTIMET_SCALAPACK_COLLECTIVES_H

#include <mpi.h>
#include <type_traits>
#include <vector>
#include "Types.h"

namespace optimet {
namespace scalapack {
namespace details {
//! True if a type scalapack knows about
template <class T>
struct is_fundamental
    : std::integral_constant<bool, std::is_same<T, int>::value or std::is_same<T, float>::value or
                                       std::is_same<T, double>::value or
                                       std::is_same<T, std::complex<float>>::value or
                                       std::is_same<T, std::complex<double>>::value> {};
}

class Context;
#define OPTIMET_MACRO(TYPE) \
  /** Broadcasts from given process to all others **/                            \
  TYPE broadcast(TYPE const&, Context const &context, t_uint row, t_uint col);   \
  Matrix<TYPE> broadcast(Matrix<TYPE> const&, Context const &context, t_uint row, t_uint col);
OPTIMET_MACRO(int);
OPTIMET_MACRO(float);
OPTIMET_MACRO(double);
OPTIMET_MACRO(std::complex<double>);
OPTIMET_MACRO(std::complex<float>);
#undef OPTIMET_MACRO

// template <class T>
// typename std::enable_if<details::is_fundamental<T>::value, T>::type
// broadcast(T const &value, Context const &context, t_uint row, t_uint col) {
//   assert(row < context.rows());
//   assert(col < context.col());
//   T result = value;
//   if(
//   OPTIMET_FC_GLOBAL_(&result, 1, registered_type(result), root, *context);
//   return result;
// }
// template <class T>
// typename std::enable_if<details::is_fundamental<T>::value, T>::type
// broadcast(Comtext const &context, t_uint row, t_uint col) {
//   assert(root < context.size());
//   assert(root != context.rank());
//   T result;
//   MPI_Bcast(&result, 1, registered_type(result), root, *context);
//   return result;
// }

// template <class T>
// Matrix<T> broadcast(Matrix<T> const &mat, Communicator const &context, t_uint root) {
//   assert(root < context.size());
//   auto const nrows = broadcast(mat.rows(), context, root);
//   auto const ncols = broadcast(mat.cols(), context, root);
//   if(context.rank() == root) {
//     MPI_Bcast(const_cast<T*>(mat.data()), nrows * ncols, Type<T>::value, root, *context);
//     return mat;
//   } else {
//     Matrix<T> result = Matrix<T>::Zero(nrows, ncols);
//     MPI_Bcast(result.data(), nrows * ncols, Type<T>::value, root, *context);
//     return result;
//   }
// }
//
// template <class MATRIX>
// typename std::enable_if<std::is_same<Matrix<typename MATRIX::Scalar>, MATRIX>::value, MATRIX>::type
// broadcast(Communicator const &context, t_uint root) {
//   assert(root < context.size());
//   assert(root != context.rank());
//   auto const nrows = broadcast<t_uint>(context, root);
//   auto const ncols = broadcast<t_uint>(context, root);
//   MATRIX result = Matrix<typename MATRIX::Scalar>::Zero(nrows, ncols);
//   MPI_Bcast(result.data(), nrows * ncols, Type<typename MATRIX::Scalar>::value, root, *context);
//   return result;
// }

} /* optime::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_MPI_COMMUNICATOR */
