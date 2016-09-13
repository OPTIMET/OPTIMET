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
//! 25, issue 4 pages 1344-1381 (2004), doi: 10.1137/s1064827501399705.
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
  Matrix<t_real> basis_rotation() const { return basis_rotation(theta(), phi(), chi()); }
  static Matrix<t_real> basis_rotation(t_real theta, t_real phi, t_real chi);

  //! \brief Simplifies access to spherical harmonics
  //! \note There is a (-1)^m factor (Condon-Shortley phase term) missing with respect to the
  //! definition used int the Gumerov paper (DOI: 10.1137/S1064827501399705).
  Complex spherical_harmonic(t_uint n, t_int m) const {
    return spherical_harmonic(n, m, theta_, phi_);
  }
  //! \brief Simplifies access to spherical harmonics
  //! \note There is a (-1)^m factor (Condon-Shortley phase term) missing with respect to the
  //! definition used int the Gumerov paper (DOI: 10.1137/S1064827501399705).
  template <class T> static std::complex<T> spherical_harmonic(t_uint n, t_int m, T theta, T phi) {
    return (m > 0 and m % 2 == 1) ? -boost::math::spherical_harmonic(n, m, theta, phi) :
                                    boost::math::spherical_harmonic(n, m, theta, phi);
  }

protected:
  //! Rotation angle in rad
  Real const theta_;
  //! Rotation angle in rad
  Real const phi_;
  //! Rotation angle in rad
  Real const chi_;

  Coefficients factors(t_uint n, t_int m, t_int mu) const;
  //! Cache which holds previously computed values
  std::map<Index, Complex> cache;
  //! Initial values
  Complex initial(t_uint n, t_int mu) const {
    return std::sqrt(4 * constant::pi / static_cast<Real>(2 * n + 1)) *
           spherical_harmonic(n, -mu, theta_, -phi_);
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

  //! \brief Performs rotation (for a single particle pair)
  //! \details The input and output consist of two-column matrices, where the first columns are the
  //! Φ coefficients, and the second columns are the Ψ coefficients. The columns should contain all
  //! coefficients up to nmax.
  template <class T0, class T1>
  void operator()(Eigen::MatrixBase<T0> const &in, Eigen::MatrixBase<T1> &out) const;
  //! \brief Performs rotation (for a single particle pair)
  //! \details The input and output consist of two-column matrices, where the first columns are the
  //! Φ coefficients, and the second columns are the Ψ coefficients. The columns should contain all
  //! coefficients up to nmax.
  template <class T0> Matrix<typename T0::Scalar> operator()(Eigen::MatrixBase<T0> const &in) const;
  //! \brief Performs adjoint/inverse rotation (for a single particle pair)
  //! \details The input and output consist of two-column matrices, where the first columns are the
  //! Φ coefficients, and the second columns are the Ψ coefficients. The columns should contain all
  //! coefficients up to nmax.
  template <class T0, class T1>
  void adjoint(Eigen::MatrixBase<T0> const &in, Eigen::MatrixBase<T1> &out) const;
  //! \brief Performs adjoint/inverse rotation (for a single particle pair)
  //! \details The input and output consist of two-column matrices, where the first columns are the
  //! Φ coefficients, and the second columns are the Ψ coefficients. The columns should contain all
  //! coefficients up to nmax.
  template <class T0> Matrix<typename T0::Scalar> adjoint(Eigen::MatrixBase<T0> const &in) const;

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

template <class T0, class T1>
void Rotation::operator()(Eigen::MatrixBase<T0> const &in, Eigen::MatrixBase<T1> &out) const {
  out.resize(in.rows(), in.cols());
  for(t_uint n(1), i(0); n < order.size(); ++n) {
    assert(in.rows() > i + order[n].rows());
    out.block(i, 0, order[n].cols(), in.rows()) =
        order[n] * in.block(i, 0, order[n].rows(), in.cols());
    i += order[n].rows();
  }
}

template <class T0>
Matrix<typename T0::Scalar> Rotation::operator()(Eigen::MatrixBase<T0> const &in) const {
  Matrix<typename T0::Scalar> out = Matrix<typename T0::Scalar>::Zero(in.rows(), in.cols());
  operator()(in, out);
  return out;
}

template <class T0, class T1>
void Rotation::adjoint(Eigen::MatrixBase<T0> const &in, Eigen::MatrixBase<T1> &out) const {
  out.resize(in.rows(), in.cols());
  for(t_uint n(1), i(0); n < order.size(); ++n) {
    assert(in.rows() > i + order[n].cols());
    out.block(i, 0, order[n].cols(), in.cols()) =
        order[n].adjoint() * in.block(i, 0, order[n].cols(), in.cols());
    i += order[n].cols();
  }
}

template <class T0>
Matrix<typename T0::Scalar> Rotation::adjoint(Eigen::MatrixBase<T0> const &in) const {
  Matrix<typename T0::Scalar> out = Matrix<typename T0::Scalar>::Zero(in.rows(), in.cols());
  adjoint(in, out);
  return out;
}
}
#endif
