#ifndef OPTIMET_ROTATION_RECURSION_H
#define OPTIMET_ROTATION_RECURSION_H

#include "Types.h"
#include "constants.h"
#include <map>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace optimet {
//! \brief Spherical harmonics projects onto rotated Spherical Harmonics
//! \details Implementation follows Nail A. Gumerov, Ramani Duraiswami, SIAM J. Sci. Comput. vol
//! 25, issue 4 pages 1344-1381 (2004), DOI: 10.1137/S1064827501399705.
class RotationCoefficients {
public:
  typedef std::tuple<t_uint, t_int, t_int> Index;
  //! Rotation angle in rad
  t_real const theta;
  //! Rotation angle in rad
  t_real const phi;
  //! Rotation angle in rad
  t_real const chi;

  RotationCoefficients(t_real const &theta, t_real const &phi, t_real const &chi)
      : theta(theta), phi(phi), chi(chi) {}

  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  //! \details from Y^m_n = \sum_\mu T_n^{\mu,n}Y^\mu_n. This operator gives T_n^{\mu, n}.
  t_complex operator()(t_uint n, t_int m, t_int mu);
  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  t_complex operator()(Index const &index) {
    return operator()(std::get<0>(index), std::get<1>(index), std::get<2>(index));
  }

  //! \brief Simplifies access to spherical harmonics
  //! \note There is a (-1)^m factor (Condon-Shortley phase term) missing with respect to the
  //! definition used int the Gumerov paper.
  t_complex spherical_harmonic(t_uint n, t_int m) const {
    return boost::math::spherical_harmonic(n, m, theta, phi);
  }

protected:
  //! Cache which holds previously computed values
  std::map<Index, t_complex> cache;
  //! Initial values
  t_complex initial(t_uint n, t_int mu) const {
    return (std::abs(mu) % 2 == 0 ? 1 : -1) *
           std::sqrt(4 * constant::pi / static_cast<t_real>(2 * n + 1)) *
           spherical_harmonic(n, -mu);
  }
  //! Applies recursion
  t_complex recursion(t_uint n, t_int m, t_int mu);
};
}
#endif
