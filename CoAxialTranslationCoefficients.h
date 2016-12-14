#ifndef COAXIAL_TRANSLATION_COEFFICIENTS_H
#define COAXIAL_TRANSLATION_COEFFICIENTS_H
#include "Types.h"
#include <array>
#include <map>
#include <type_traits>
#include <iostream>
#include <vector>

#include "Spherical.h"

namespace optimet {
class CachedCoAxialRecurrence {
public:
  class Functor {
  public:
    //! Creates from coefficients that are moved here
    Functor(t_int N, std::vector<t_complex> &&coeffs) : N(N), coefficients(std::move(coeffs)) {}
    //! Applies direct functor
    template <class T0, class T1>
    typename std::enable_if<std::is_same<typename T0::Scalar, t_complex>::value>::type
    operator()(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out) const;
    //! Applies direct functor
    template <class T>
    typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
    operator()(Eigen::MatrixBase<T> const &input) const;
    //! Applies transpose functor
    template <class T0, class T1>
    typename std::enable_if<std::is_same<typename T0::Scalar, t_complex>::value>::type
    transpose(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out) const;
    //! Applies transpose functor
    template <class T>
    typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
    transpose(Eigen::MatrixBase<T> const &input) const;

  private:
    t_int N;
    std::vector<t_complex> coefficients;
  };
  //! Inner floating point with higher precision
  typedef long double Real;
  //! Inner complex floating point with higher precision
  typedef std::complex<Real> Complex;
  //! Indices tuple
  typedef std::array<t_int, 3> t_indices;

  CachedCoAxialRecurrence(t_real distance, t_complex waveK, bool regular = true)
      : distance(distance), waveK(waveK), regular(regular) {}

  //! \brief Returns coaxial translation coefficients
  //! \details n, l and m correspond to the same variables in Gumerov (2002),
  //! s = m by definition.
  t_complex operator()(t_int n, t_int m, t_int l) { return static_cast<t_complex>(coeff(n, m, l)); }

  bool is_regular() const { return regular; }

  //! \brief Applies recurrence to input vector/matrix
  //! \details Each input column consists of (n, m) elements arranged in descending order (1, -1),
  //! (1, 0), (1, 1), (2, -2), ... (nmax, nmax). nmax is determined from the number of rows.
  template <class T0, class T1>
  typename std::enable_if<std::is_same<typename T0::Scalar, t_complex>::value>::type
  operator()(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out);

  //! \brief Applies recurrence to input vector/matrix
  //! \details Each input column consists of (n, m) elements arranged in descending order (1, -1),
  //! (1, 0), (1, 1), (2, -2), ... (nmax, nmax). nmax is determined from the number of rows.
  template <class T>
  typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
  operator()(Eigen::MatrixBase<T> const &input);

  // \brief creates a functor that can be used to apply the matrix
  // \details The functor takes an input and output vector (of radiating or non-radiating
  // coeffiecients). The ouput vector contains the result of
  // applying the coaxial translation to the input vector. The functor is specialized for a
  // specific number of harmonics/input vector size. It is an error to call it on a vector with a
  // different size.
  Functor functor(t_int n);

protected:
  //! Distance that the solution is to be translated by
  Real const distance;
  //! Wavenumber of the incident wave
  Complex const waveK;
  //! Whether this is for regular or irregular coeffs
  bool const regular;
  //! Caches known coefficients
  std::map<t_indices, Complex> cache;

  //! Switches between recurrence relationships
  Complex recurrence(t_int n, t_int m, t_int l);
  Complex initial(t_int l);
  Complex sectorial_recurrence(t_int n, t_int m, t_int l);
  Complex zonal_recurrence(t_int n, t_int l);
  Complex offdiagonal_recurrence(t_int n, t_int m, t_int l);
  Complex coeff(t_int n, t_int m, t_int l);
};

template <class T0, class T1>
typename std::enable_if<std::is_same<typename T0::Scalar, t_complex>::value>::type
CachedCoAxialRecurrence::
operator()(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out) {
  auto const nr = input.rows();
  auto const with_n0 = std::abs(std::sqrt(nr) - std::lround(std::sqrt(nr))) <
                       std::abs(std::sqrt(nr + 1) - std::lround(std::sqrt(nr + 1)));
  t_int const N = std::lround(std::sqrt(with_n0 ? nr : nr + 1)) - 1;
  auto const min_n = with_n0 ? 0 : 1;
  assert((with_n0 and (N + 1) * (N + 1) == input.rows()) or N * (N + 2) == input.rows());
  const_cast<Eigen::MatrixBase<T1> &>(out).resize(input.rows(), input.cols());
  const_cast<Eigen::MatrixBase<T1> &>(out).fill(0);
  assert(N >= min_n);
  auto const index = with_n0 ? [](t_int n, t_int m) { return n * (n + 1) + m; } :
                               [](t_int n, t_int m) { return n * (n + 1) + m - 1; };
  assert(index(min_n, -min_n) == 0);
  assert(index(N, N) + 1 == input.rows());
  for(auto n = min_n, i = 0; n <= N; ++n)
    for(auto m = -n; m <= n; ++m, ++i)
      for(auto l = std::max(std::abs(m), min_n); l <= N; ++l)
        const_cast<Eigen::MatrixBase<T1> &>(out).row(index(l, m)) += operator()(n, m, l) *
                                                                     input.row(i);
}

template <class T>
typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
CachedCoAxialRecurrence::operator()(Eigen::MatrixBase<T> const &input) {
  typedef typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>,
                                    Matrix<t_complex>>::type Out;
  Out out(input.rows(), input.cols());
  operator()(input, out);
  return out;
}

template <class T0, class T1>
typename std::enable_if<std::is_same<typename T0::Scalar, t_complex>::value>::type
CachedCoAxialRecurrence::Functor::
operator()(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out) const {
  auto const nr = input.rows();
  auto const with_n0 = std::abs(std::sqrt(nr) - std::lround(std::sqrt(nr))) <
                       std::abs(std::sqrt(nr + 1) - std::lround(std::sqrt(nr + 1)));
  t_int const N = std::lround(std::sqrt(with_n0 ? nr : nr + 1)) - 1;
  if(N > this->N)
    throw std::out_of_range("Input matrix too large");
  assert((with_n0 and (N + 1) * (N + 1) == input.rows()) or N * (N + 2) == input.rows());
  int const min_n = with_n0 ? 0 : 1;
  auto const index = with_n0 ? [](t_int n, t_int m) { return n * (n + 1) + m; } :
                               [](t_int n, t_int m) { return n * (n + 1) + m - 1; };
  assert(index(min_n, -min_n) == 0);
  assert(index(N, N) + 1 == input.rows());
  const_cast<Eigen::MatrixBase<T1> &>(out).resize(input.rows(), input.cols());
  const_cast<Eigen::MatrixBase<T1> &>(out).fill(0);
  auto i_coeff = with_n0 ? coefficients.begin() : (coefficients.begin() + N + 1);
  for(auto n = min_n, i = 0; n <= N; ++n)
    for(auto m = -n; m <= n; ++m, ++i)
      for(auto l = std::abs(m); l <= N; ++l, ++i_coeff) {
        assert(i_coeff != coefficients.end());
        auto const i_out = index(l, m);
        assert(i < input.rows());
        if(i_out >= 0 and i_out < out.rows())
          const_cast<Eigen::MatrixBase<T1> &>(out).row(i_out) += (*i_coeff) * input.row(i);
      }
}

template <class T0, class T1>
typename std::enable_if<std::is_same<typename T0::Scalar, t_complex>::value>::type
CachedCoAxialRecurrence::Functor::transpose(Eigen::MatrixBase<T0> const &input,
                                            Eigen::MatrixBase<T1> const &out) const {
  auto const nr = input.rows();
  auto const with_n0 = std::abs(std::sqrt(nr) - std::lround(std::sqrt(nr))) <
                       std::abs(std::sqrt(nr + 1) - std::lround(std::sqrt(nr + 1)));
  t_int const N = std::lround(std::sqrt(with_n0 ? nr : nr + 1)) - 1;
  if(N > this->N)
    throw std::out_of_range("Input matrix too large");
  assert((with_n0 and (N + 1) * (N + 1) == input.rows()) or N * (N + 2) == input.rows());
  int const min_n = with_n0 ? 0 : 1;
  auto const index = with_n0 ? [](t_int n, t_int m) { return n * (n + 1) + m; } :
                               [](t_int n, t_int m) { return n * (n + 1) + m - 1; };
  assert(index(min_n, -min_n) == 0);
  assert(index(N, N) + 1 == input.rows());
  const_cast<Eigen::MatrixBase<T1> &>(out).resize(input.rows(), input.cols());
  const_cast<Eigen::MatrixBase<T1> &>(out).fill(0);
  auto i_coeff = with_n0 ? coefficients.begin() : (coefficients.begin() + N + 1);
  for(auto n = min_n, i = 0; n <= N; ++n)
    for(auto m = -n; m <= n; ++m, ++i)
      for(auto l = std::abs(m); l <= N; ++l, ++i_coeff) {
        assert(i_coeff != coefficients.end());
        assert(i < out.rows());
        auto const i_in = index(l, m);
        if(i_in >= 0 and i_in < input.rows())
          const_cast<Eigen::MatrixBase<T1> &>(out).row(i) += (*i_coeff) * input.row(i_in);
      }
}

template <class T>
typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
CachedCoAxialRecurrence::Functor::operator()(Eigen::MatrixBase<T> const &input) const {
  typedef typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>,
                                    Matrix<t_complex>>::type Out;
  Out out(input.rows(), input.cols());
  operator()(input, out);
  return out;
}

template <class T>
typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
CachedCoAxialRecurrence::Functor::transpose(Eigen::MatrixBase<T> const &input) const {
  typedef typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>,
          Matrix<t_complex>>::type Out;
  Out out(input.rows(), input.cols());
  transpose(input, out);
  return out;
}
}
#endif
