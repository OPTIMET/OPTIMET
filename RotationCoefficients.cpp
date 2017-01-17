#include "RotationCoefficients.h"
#include <Coefficients.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iterator>

namespace optimet {

RotationCoefficients::Complex RotationCoefficients::with_caching(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return static_cast<Real>(0);
  if(m < 0)
    return std::conj(with_caching(n, -m, -mu));
  if(n == 0)
    return 1;

  auto const prior = cache.find(std::make_tuple(n, m, mu));
  if(prior != cache.end())
    return prior->second;

  auto const value = recursion(n, m, mu);
  cache[std::make_tuple(n, m, mu)] = value;
  return value;
}

RotationCoefficients::Complex RotationCoefficients::recursion(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return static_cast<Real>(0);
  if(m < 0)
    return std::conj(with_caching(n, -m, -mu));

  if(m == 0)
    return initial(n, mu);

  using coefficient::a;
  using coefficient::b;

  auto const factor = std::exp(Complex(0, chi_)) / b<Real>(n + 1, m - 1);
  auto const half_factor = static_cast<Real>(0.5) * factor;
  auto const c0 =
      half_factor * b<Real>(n + 1, -mu - 1) * std::exp(Complex(0, phi_)) * (1 - std::cos(theta_));
  auto const c1 =
      -half_factor * b<Real>(n + 1, mu - 1) * std::exp(Complex(0, -phi_)) * (1 + std::cos(theta_));
  auto const c2 = -factor * a<Real>(n, mu) * std::sin(theta_);

  return c0 * with_caching(n + 1, m - 1, mu + 1) + c1 * with_caching(n + 1, m - 1, mu - 1) +
         c2 * with_caching(n + 1, m - 1, mu);
}

Matrix<t_complex> RotationCoefficients::matrix(t_uint n) {
  Matrix<t_complex> result = Matrix<t_complex>::Zero(2 * n + 1, 2 * n + 1);
  for(t_uint i(0); i < 2 * n + 1; ++i) {
    t_int const m = static_cast<t_int>(i) - static_cast<t_int>(n);
    for(t_uint j(0); j < 2 * n + 1; ++j) {
      t_int const mu = static_cast<t_int>(j) - static_cast<t_int>(n);
      result(i, j) = operator()(n, m, mu);
    }
  }
  return result;
}

Rotation::Rotation(t_real const &theta, t_real const &phi, t_real const &chi, t_uint nmax)
    : theta_(theta), phi_(phi), chi_(chi), nmax_(nmax) {
  order.reserve(nmax + 1);
  RotationCoefficients coeffs(theta, phi, chi);
  for(t_uint i(0); i <= nmax; ++i)
    order.push_back(coeffs.matrix(i));
}

Eigen::Matrix<t_real, 3, 3>
    RotationCoefficients::basis_rotation(Eigen::Matrix<t_real, 3, 1> const &axis) {
  if(axis.stableNorm() < 1e-8)
    throw std::runtime_error("Input axis is zero");

  Eigen::Matrix<t_real, 3, 1> const zhat = axis.normalized();
  if(zhat.isApprox(zhat.Unit(2)))
    return Eigen::Matrix<t_real, 3, 3>::Identity();
  if(zhat.isApprox(-zhat.Unit(2))) {
    Eigen::Matrix<t_real, 3, 3> result;
    result.fill(0);
    result.diagonal<0>() = Eigen::Matrix<t_real, 3, 1>(1, -1, -1);
    return result;
  }

  auto const xhat = zhat.Unit(2).cross(zhat).normalized().eval();
  auto const yhat = zhat.cross(xhat).normalized().eval();
  Eigen::Matrix<t_real, 3, 3> result;
  result.row(0) = xhat;
  result.row(1) = yhat;
  result.row(2) = zhat;
  return result;
}

Eigen::Matrix<t_real, 3, 3>
RotationCoefficients::basis_rotation(t_real theta, t_real phi, t_real chi) {
  using std::sin;
  using std::cos;
  Eigen::Matrix<t_real, 3, 3> result;
  result.col(0) << -sin(phi) * sin(chi) - cos(theta) * cos(phi) * cos(chi),
      cos(phi) * sin(chi) - cos(theta) * sin(phi) * cos(chi), sin(theta) * cos(chi);
  result.col(1) << sin(phi) * cos(chi) - cos(theta) * cos(phi) * sin(chi),
      -(cos(phi) * cos(chi) + cos(theta) * sin(phi) * sin(chi)), sin(theta) * sin(chi);
  result.col(2) << sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta);
  return result;
}

std::tuple<t_real, t_real, t_real>
    RotationCoefficients::rotation_angles(Eigen::Matrix<t_real, 3, 3> const &matrix) {
  auto const is_zaxis = std::abs(std::abs(matrix(2, 2)) - 1) < 1e-12;
  auto const theta = std::acos(matrix(2, 2));
  auto const phi = is_zaxis ? 0 : std::atan2(matrix(1, 2), matrix(0, 2));
  auto const chi = is_zaxis ?
                       std::atan2(matrix(1, 0), matrix(2, 2) > 0 ? -matrix(0, 0) : matrix(0, 0)) :
                       std::atan2(matrix(2, 1), matrix(2, 0));
  if(!RotationCoefficients::basis_rotation(theta, phi, chi).isApprox(matrix))
    throw std::runtime_error("Could not recover rotation matrix");
  return std::tuple<t_real, t_real, t_real>(theta, phi, chi);
}
}
