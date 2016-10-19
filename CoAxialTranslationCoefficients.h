#ifndef COAXIAL_TRANSLATION_COEFFICIENTS_H
#include "Types.h"
#include <array>
#include <iostream>
#include <map>
#include <type_traits>

#include "Spherical.h"

namespace optimet {
class CachedCoAxialRecurrence {
public:
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
  operator()(Eigen::MatrixBase<T0> &out, Eigen::MatrixBase<T1> const &input);

  //! \brief Applies recurrence to input vector/matrix
  //! \details Each input column consists of (n, m) elements arranged in descending order (1, -1),
  //! (1, 0), (1, 1), (2, -2), ... (nmax, nmax). nmax is determined from the number of rows.
  template <class T>
  typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>, Matrix<t_complex>>::type
  operator()(Eigen::MatrixBase<T> const &input) {
    typedef typename std::conditional<T::ColsAtCompileTime == 1, Vector<t_complex>,
                                      Matrix<t_complex>>::type Out;
    Out out(input.rows(), input.cols());
    operator()(out, input);
    return out;
  }

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
operator()(Eigen::MatrixBase<T0> &out, Eigen::MatrixBase<T1> const &input) {
  out.resize(input.rows(), input.cols());
  out.fill(0);
  t_int const nmax = std::lround(std::sqrt(input.rows()) - 1.0);
  assert((nmax + 1) * (nmax + 1) == input.rows());
  assert(nmax >= 0);
  auto const index = [](t_int n, t_int m) { return std::abs(m) > n ? 0 : n * n + n + m; };
  for(auto n = 0, i = 0; n <= nmax; ++n)
    for(auto m = -n; m <= n; ++m, ++i)
      for(auto l = std::abs(m); l <= nmax; ++l) {
        assert(i == index(n, m));
        out.row(index(l, m)) += operator()(n, m, l) * input.row(index(n, m));
      }
}
}
#endif
