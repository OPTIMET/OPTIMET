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
  //! Inner floating point with higher precision
  typedef long double Real;
  //! Inner complex floating point with higher precision
  typedef std::complex<Real> Complex;

public:
  typedef std::tuple<t_uint, t_int, t_int> Index;
  typedef Eigen::Matrix<Complex, 3, 1> Coefficients;

  RotationCoefficients(t_real const &theta, t_real const &phi, t_real const &chi)
      : theta_(static_cast<Real>(theta)), phi_(static_cast<Real>(phi)),
        chi_(static_cast<Real>(chi)) {}

  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  //! \details from Y^m_n = \sum_\mu T_n^{\mu,n}Y^\mu_n. This operator gives T_n^{\mu, n}.
  t_complex operator()(t_uint n, t_int m, t_int mu) {
    return static_cast<t_complex>(with_caching(n, m, mu));
  }
  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  //! \brief Spherical Harmonic Y^m_n projected onto Y^\mu_n
  t_complex operator()(Index const &index) {
    return operator()(std::get<0>(index), std::get<1>(index), std::get<2>(index));
  }

  t_real theta() const { return static_cast<t_real>(theta_); }
  t_real phi() const { return static_cast<t_real>(phi_); }
  t_real chi() const { return static_cast<t_real>(chi_); }

  //! Matrix for spherical harmonics of given degree
  Matrix<t_complex> matrix(t_uint n);

  //! Basis rotation from a to ^a
  Matrix<t_real> basis_rotation() const {
    return basis_rotation(theta(), phi(), chi());
  }
  static Matrix<t_real> basis_rotation(t_real theta, t_real phi, t_real chi);

protected:
  //! Rotation angle in rad
  Real const theta_;
  //! Rotation angle in rad
  Real const phi_;
  //! Rotation angle in rad
  Real const chi_;

  //! \brief Simplifies access to spherical harmonics
  //! \note There is a (-1)^m factor (Condon-Shortley phase term) missing with respect to the
  //! definition used int the Gumerov paper.
  Complex spherical_harmonic(t_uint n, t_int m) const {
    return (m > 0 and m % 2 == 1) ? -boost::math::spherical_harmonic(n, m, theta_, phi_) :
                                    boost::math::spherical_harmonic(n, m, theta_, phi_);
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

//! \brief Rotation by (phi, psi, chi) for orders up to nmax
class Rotation {
   public:
     Rotation(t_real const &theta, t_real const &phi, t_real const &chi, t_uint nmax);

     t_real theta() const { return theta_; }
     t_real phi() const { return phi_; }
     t_real chi() const { return chi_; }
     t_uint nmax() const { return nmax_; }

     //! creates a rotation matrix for the given input
     Matrix<t_complex> rotation_matrix(t_real n) {
       return rotation_matrix(RotationCoefficients(theta(), phi(), chi()), n);
     }
     //! creates a rotation matrix for the given input
     Matrix<t_complex> rotation_matrix(RotationCoefficients &coeffs, t_real n) {
       return coeffs.matrix(n);
     }
     //! creates a rotation matrix for the given input
     Matrix<t_complex> rotation_matrix(RotationCoefficients &&coeffs, t_real n) {
       return coeffs.matrix(n);
     }

   protected:
     //! Rotation angle in rad
     t_real const theta_;
     //! Rotation angle in rad
     t_real const phi_;
     //! Rotation angle in rad
     t_real const chi_;
     //! Maximum degree of the spherical harmonics
     t_uint const nmax_;
     //! Matrices for each spherical harmonic up to given order
     std::vector<Matrix<t_complex>> order;
};
}
#endif
