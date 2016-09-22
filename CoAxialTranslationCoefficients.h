#ifndef COAXIAL_TRANSLATION_COEFFICIENTS_H
#include "Types.h"
#include <array>
#include <iostream>
#include <map>

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
  Complex operator()(t_int n, t_int m, t_int l);

  bool is_regular() const { return regular; }

  //! \brief Applies recurrence to input vector/matrix
  //! \details Each input column consists of (n, m) elements arranged in descending order (1, -1),
  //! (1, 0), (1, 1), (2, -2), ... (nmax, nmax). nmax is determined from the number of rows.
  template <class T0, class T1>
  void operator()(Eigen::MatrixBase<T0> &out, Eigen::MatrixBase<T1> const &input);

  //! \brief Applies recurrence to input vector/matrix
  template <class T> Matrix<typename T::Scalar> operator()(Eigen::MatrixBase<T> const &input) {
    Matrix<typename T::Scalar> out(input.rows(), input.cols());
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
};

template <class T0, class T1>
void CachedCoAxialRecurrence::
operator()(Eigen::MatrixBase<T0> &out, Eigen::MatrixBase<T1> const &input) {

  out.resize(input.rows(), input.cols());
  t_int const nmax = std::lround(std::sqrt(input.rows() + 1) - 1.0);
  assert(nmax * (nmax + 2) == input.rows());
  assert(nmax > 0);
  auto const index = [](t_int n, t_int m) {
    // for |m| > n, return something valid, e.g. 0.
    // (n - 1) * (n + 1) --> same as nmax, but for n - 1 --> n * n - 1
    // n + m => index m in [-n, n] becomes index in [0, 2n]
    return std::abs(m) > n ? 0 : n * n - 1 + n + m;
  };

  std::cout << "rows: " << input.rows() << ", " << out.rows() << "\n";
  for(auto n = 1; n <= nmax; ++n)
    for(auto m = -n; m <= n; ++m) {
      auto current_row = out.row(index(n, m));
      current_row.fill(0);
      for(auto l = std::max(1, std::abs(m)); l <= nmax; ++l)
        current_row.array() += operator()(n, m, l) * input.row(index(l, m)).array();
    }
}
}
#endif
