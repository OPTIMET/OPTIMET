#ifndef OPTIMET_ROTATION_RECURSION_H
#define OPTIMET_ROTATION_RECURSION_H

#include "Types.h"
#include "constants.h"
#include <map>
#include <vector>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace optimet {
//! \brief Spherical harmonics projects onto rotated Spherical Harmonics
//! \details Implementation follows Nail A. Gumerov, Ramani Duraiswami, SIAM J. Sci. Comput. vol
//! 25, issue 4 pages 1344-1381 (2004), DOI: 10.1137/S1064827501399705.
class RotationCoefficients {
public:
  typedef t_complex Complex;
  typedef Complex::value_type Real;

  typedef std::tuple<t_uint, t_int, t_int> Index;
  typedef Eigen::Matrix<Complex, 3, 1> Coefficients;
  //! Rotation angle in rad
  Real const theta;
  //! Rotation angle in rad
  Real const phi;
  //! Rotation angle in rad
  Real const chi;

  RotationCoefficients(Real const &theta, Real const &phi, Real const &chi)
      : theta(theta), phi(phi), chi(chi) {}

  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  //! \details from Y^m_n = \sum_\mu T_n^{\mu,n}Y^\mu_n. This operator gives T_n^{\mu, n}.
  Complex operator()(t_uint n, t_int m, t_int mu) {
    return with_caching(n, m, mu);
  }
  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  Complex operator()(Index const &index) {
    return operator()(std::get<0>(index), std::get<1>(index), std::get<2>(index));
  }

protected:
  //! \brief Simplifies access to spherical harmonics
  //! \note There is a (-1)^m factor (Condon-Shortley phase term) missing with respect to the
  //! definition used int the Gumerov paper.
  Complex spherical_harmonic(t_uint n, t_int m) const {
    return (m > 0 and m % 2 == 1) ? -boost::math::spherical_harmonic(n, m, theta, phi) :
                                    boost::math::spherical_harmonic(n, m, theta, phi);
  }

  static Real a(t_uint n, t_int m);
  static Real b(t_uint n, t_int m);
  Coefficients factors(t_uint n, t_int m, t_int mu) const;
  //! Cache which holds previously computed values
  std::map<Index, Complex> cache;
  //! Initial values
  Complex initial(t_uint n, t_int mu) const {
    return std::sqrt(4 * constant::pi / static_cast<Real>(2 * n + 1)) * spherical_harmonic(n, -mu);
  }
  //! Applies recursion
  Complex recursion(t_uint n, t_int m, t_int mu);
  //! Applies caching to recursion
  Complex with_caching(t_uint n, t_int m, t_int mu);
};
}
#endif
