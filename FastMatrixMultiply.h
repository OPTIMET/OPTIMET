#ifndef FAST_MATRIX_MULTIPLY_H

#include "Bessel.h"
#include "CoAxialTranslationCoefficients.h"
#include "RotationCoaxialDecomposition.h"
#include "RotationCoefficients.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Types.h"
#include <utility>
#include <vector>

namespace optimet {
//! Electric and Magnetic field coefficients, from scalar potentials of a single particle
template <class T0, class T1>
void field_coefficients(t_real wavenumber, Eigen::MatrixBase<T0> const &coeffs,
                        Eigen::MatrixBase<T1> const &out);

//! Computes fields at a given position
class Fields {
public:
  //! Position in cartesian coordinates
  typedef Eigen::Matrix<t_real, 3, 1> Position;
  //! First column holds E, second H
  typedef Eigen::Matrix<t_complex, 3, 2> Result;

  template <class T>
  Fields(std::vector<Scatterer> const &objects, t_real wavenumber, ElectroMagnetic const &bground,
         Eigen::DenseBase<T> const &coeffs)
      : outside(coeffs.size(), 6), inside(coeffs.size(), 6), objects(objects),
        wavenumber(wavenumber) {
    auto const coutside = convertIndirect(coeffs, wavenumber * constant::c, bground, objects);
    auto const cinside = convertInternal(outside, wavenumber * constant::c, bground, objects);

    t_int row(0);
    for(auto const &object : objects) {
      auto const N = 2 * object.nMax * (object.nMax + 2);
      field_coefficients(wavenumber, coutside.segment(row, N), outside.middleRows(row, N));
      field_coefficients(wavenumber, cinside.segment(row, N), inside.middleRows(row, N));
      row += N;
    }
  }

  //! Computes field at given position
  Result operator()(Position const &position) const;

protected:
  Matrix<t_complex> outside;
  Matrix<t_complex> inside;
  std::vector<Scatterer> objects;
  t_real wavenumber;
  Result compute_outside(Position const &position) const;
  template <class T>
  Result
  single_sphere(Position const &position, t_int nMax, Eigen::MatrixBase<T> const &coeffs) const;
};

//! Returns function to compute the fields at a given position
std::function<Eigen::Matrix<t_complex, 3, 2>(Eigen::Matrix<t_real, 3, 1> const &)>
fields(t_real wavenumber, Vector<t_complex> const &coeffs, std::vector<Scatterer> const &objects);

//! \brief Fast multiplication of effective incident field by transfer matrix
//! \details The multiplication occurs for a set of scattering particles receiving the EM
//! radiation
//! to a set of locations where the fields are checked. Each incident effective field is
//! composed of coefficients for both Φ and Ψ potentials (in R basis).
class FastMatrixMultiply {
public:
  //! Range of incident/scattering objects
  typedef std::pair<t_int, t_int> Range;

  //! Creates the fast matrix multiply object
  //! \param[in] em_background: Electro-magnetic properties of the background material
  //! \param[in] wavenumber: Angular wave-number of the incident plane-wave
  //! \param[in] scatterers: all spherical scatterers in the problem
  //! \param[in] incident: Indices defining the scatterers this instance will take care of.
  //!     In practice, the input vector to the multiplication will contain the coefficients of the
  //!     effective incident field expanded on a spherical basis set at the location of the
  //!     scatterers in this range.
  //! \param[in] translate: Indices defining the location for which to expand the field on output.
  //!     In practice, the output vector of the multiplication will contain the coefficients of
  //!     the
  //!     spherical basis set used to expand the field at the location of the scatterers in this
  //!     range.
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Range incident, Range translate)
      : em_background_(em_background), wavenumber_(wavenumber), scatterers_(scatterers),
        global_indices_(compute_indices(scatterers_)), incident_range_(incident),
        translate_range_(translate),
        rotations_(compute_rotations(scatterers_, incident_range_, translate_range_)),
        mie_coefficients_(
            compute_mie_coefficients(em_background, wavenumber, scatterers_, incident_range_)),
        coaxial_translations_(compute_coaxial_translations(wavenumber_, scatterers_,
                                                           incident_range_, translate_range_)) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers, Range incident,
                     Range translate)
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers, incident, translate) {}
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers)
      : FastMatrixMultiply(em_background, wavenumber, scatterers, {0, scatterers.size()},
                           {0, scatterers.size()}) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers)
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers) {}

  //! Computes index of each particle i in global input vector
  static std::vector<t_uint> compute_indices(std::vector<Scatterer> const &scatterers);
  //! Computes rotations between relevant pairs of particles
  static std::vector<Rotation>
  compute_rotations(std::vector<Scatterer> const &scatterers, Range incident, Range translate);
  //! Computes co-axial translations between relevant pairs of particles
  static std::vector<CachedCoAxialRecurrence::Functor>
  compute_coaxial_translations(t_complex wavenumber_, std::vector<Scatterer> const &scatterers,
                               Range incident, Range translate);
  //! Computes mie coefficient for each particles
  static Vector<t_complex>
  compute_mie_coefficients(ElectroMagnetic const &background, t_real wavenumber,
                           std::vector<Scatterer> const &scatterers, Range incident);
  //! Total size of the problem
  t_uint size() const { return global_indices_.back(); }

  //! \brief Applies fast matrix multiplication to effective incident field
  void operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator()(Vector<t_complex> const &in) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator*(Vector<t_complex> const &in) const { return operator()(in); }

  //! \brief Adds co-axial rotation/translation for a given particle pair
  //! \details the coefficients from input are projected from iScatt to oScatt using the given
  //! rotation translation functors. The result are added to the out-going coefficients.
  template <class T0, class T1>
  void remove_translation(Eigen::MatrixBase<T0> const &input, Scatterer const &iScatt,
                          Scatterer const &oScatt, Rotation const &rotation,
                          CachedCoAxialRecurrence::Functor const &translation,
                          Eigen::MatrixBase<T1> const &out) const;

  //! Adds co-axial rotation/translation for a given particle pair
  template <class T0, class T1, class T2>
  void
  remove_translation(Eigen::MatrixBase<T0> const &input, Scatterer const &iScatt,
                     Scatterer const &oScatt, Rotation const &rotation,
                     CachedCoAxialRecurrence::Functor const &translation,
                     Eigen::MatrixBase<T1> const &out, Eigen::MatrixBase<T2> const &work) const;

  //! Apply translation to each particle pair
  void translation(Vector<t_complex> const &in, Vector<t_complex> &out) const;

protected:
  static int const nplus = 1;
  //! The properties of the background
  ElectroMagnetic em_background_;
  //! Wavenumber of the incident wave
  t_real wavenumber_;
  //! Scatterers for which to compute matrix product
  std::vector<Scatterer> const scatterers_;
  //! Starting index of particle i
  std::vector<t_uint> const global_indices_;
  //! \brief Index of first owned sphere
  //! \details Owned spheres are those for which this object will apply Matrix-Vector
  //! multiplications
  Range incident_range_;
  //! \brief Index of first owned sphere
  //! \details Owned spheres are those for which this object will apply Matrix-Vector
  //! multiplications
  Range translate_range_;
  //! Rotations for owned objects
  std::vector<Rotation> const rotations_;
  //! Mie coefficients
  Vector<t_complex> const mie_coefficients_;
  //! Co-axial translations
  std::vector<CachedCoAxialRecurrence::Functor> coaxial_translations_;

  //! Number of basis function for given nmax
  static constexpr t_int nfunctions(t_int nmax) { return nmax * (nmax + 2); }
};

template <class T0, class T1>
void FastMatrixMultiply::remove_translation(Eigen::MatrixBase<T0> const &input,
                                            Scatterer const &iScatt, Scatterer const &oScatt,
                                            Rotation const &rotation,
                                            CachedCoAxialRecurrence::Functor const &translation,
                                            Eigen::MatrixBase<T1> const &out) const {
  if(input.cols() != 2)
    throw std::runtime_error("input should have two columns");
  if(out.cols() != 2)
    throw std::runtime_error("output should have two columns");

  auto const nmax = std::max(iScatt.nMax, oScatt.nMax) + nplus;
  Eigen::Matrix<t_complex, Eigen::Dynamic, 4> work(nfunctions(nmax) + 1, 4);
  FastMatrixMultiply::remove_translation(input, iScatt, oScatt, rotation, translation, out, work);
}

template <class T0, class T1, class T2>
void FastMatrixMultiply::remove_translation(Eigen::MatrixBase<T0> const &input,
                                            Scatterer const &iScatt, Scatterer const &oScatt,
                                            Rotation const &rotation,
                                            CachedCoAxialRecurrence::Functor const &translation,
                                            Eigen::MatrixBase<T1> const &out,
                                            Eigen::MatrixBase<T2> const &work) const {
  auto const in_rows = nfunctions(iScatt.nMax);
  auto const out_rows = nfunctions(oScatt.nMax);

  auto const max_rows = nfunctions(std::max(iScatt.nMax, oScatt.nMax) + nplus);
  if(work.rows() <= max_rows or work.cols() != 4)
    const_cast<Eigen::MatrixBase<T2> &>(work).resize(max_rows + 1, 4);
  const_cast<Eigen::MatrixBase<T2> &>(work).fill(0);

  assert((iScatt.vR.toEigenCartesian() - oScatt.vR.toEigenCartesian()).stableNorm() >=
         iScatt.radius + oScatt.radius);

  // First apply rotation - without n=0 term
  auto rotated = const_cast<Eigen::MatrixBase<T2> &>(work).leftCols(2);
  rotation.transpose(input, rotated.middleRows(1, in_rows));

  // Then add normalization coming from ???
  // But appears when comparing to original calculation
  // We add it second since: (i) the normalization is same for the same n, (ii) it avoids a copy
  for(t_int j(1), n(1); n <= iScatt.nMax; j += 2 * n + 1, ++n) {
    rotated.col(0).segment(j, 2 * n + 1) /= std::sqrt((n * (n + 1)) / 2);
    rotated.col(1).segment(j, 2 * n + 1) /= -std::sqrt((n * (n + 1)) / 2);
  }

  // Then perform co-axial translation - this may create n=0 term
  auto ztrans = const_cast<Eigen::MatrixBase<T2> &>(work).rightCols(2);
  translation(rotated, ztrans);

  // Then apply field-coaxial-tranlation transform thing - n=0 term may be used to create n=1 term.
  // n=0 term itself becomes zero (thereby choosing a gauge, apparently)
  auto const tz = (oScatt.vR.toEigenCartesian() - iScatt.vR.toEigenCartesian()).stableNorm();
  auto &field_translated = rotated;
  rotation_coaxial_decomposition(wavenumber_, tz, ztrans, field_translated);

  // Rotate back - remove n=0 term since it is zero
  auto unrotated = const_cast<Eigen::MatrixBase<T2> &>(work).block(1, 2, out_rows, 2);
  rotation.conjugate(field_translated.middleRows(1, out_rows), unrotated);

  // Add normalization coming from ???
  // But appears when comparing to original calculation
  for(t_int j(0), n(1); n <= oScatt.nMax; j += 2 * n + 1, ++n) {
    unrotated.col(0).segment(j, 2 * n + 1) *= std::sqrt((n * (n + 1)) / 2);
    unrotated.col(1).segment(j, 2 * n + 1) *= -std::sqrt((n * (n + 1)) / 2);
  }

  // Finally, add back into output vector
  const_cast<Eigen::MatrixBase<T1> &>(out) -= unrotated;
}

template <class T0, class T1>
void field_coefficients(t_real wavenumber, Eigen::MatrixBase<T0> const &coeffs,
                        Eigen::MatrixBase<T1> const &out) {
  using coefficient::a;
  using coefficient::b;
  using coefficient::c;

  assert(coeffs.size() % 2 == 0);
  t_int const N = std::lround(std::sqrt(coeffs.size() / 2 + 1)) - 1;
  assert(N * (N + 2) == coeffs.size());
  auto const index = [](t_int n, t_int m) { return std::abs(m) > n ? 0 : n * (n + 1) + m - 1; };
  assert(index(0, 0) == 0);
  assert(index(N, N) + 1 == coeffs.rows());

  const_cast<Eigen::MatrixBase<T1> &>(out).resize(coeffs.rows(), 6);
  auto E = const_cast<Eigen::MatrixBase<T1> &>(out).leftCols(3);
  auto H = const_cast<Eigen::MatrixBase<T1> &>(out).rightCols(3);

#define OPTIMET_X(PHI, PSI)                                                                        \
  t_complex(0, 0.5) * (c(n, m - 1) * PHI(n, m - 1) + c(n, m) * PHI(n, m)) -                        \
      wavenumber / 2 * ((n + 2) * b(n + 1, m - 1) * PSI(n + 1, m - 1) +                            \
                        (n - 1) * b(n, -m) * PSI(n - 1, m - 1) +                                   \
                        (n + 2) * b(n + 1, -m - 1) * PSI(n + 1, m + 1) +                           \
                        (n - 1) * b(n, m) * PSI(n - 1, m + 1));
#define OPTIMET_Y(PHI, PSI)                                                                        \
  -0.5 * (-c(n, m - 1) * PHI(n, m - 1) + c(n, m) * PHI(n, m)) +                                    \
      wavenumber *t_complex(0, 0.5) * ((n + 2) * b(n + 1, m - 1) * PSI(n + 1, m - 1) +             \
                                       (n - 1) * b(n, -m) * PSI(n - 1, m - 1) -                    \
                                       (n + 2) * b(n + 1, -m - 1) * PSI(n + 1, m + 1) -            \
                                       (n - 1) * b(n, m) * PSI(n - 1, m + 1));
#define OPTIMET_Z(PHI, PSI)                                                                        \
  t_complex(0, -m) * PHI(n, m) +                                                                   \
      wavenumber *((n + 2) * a(n, m) * PSI(n + 1, m) + (n - 1) * a(n - 1, m) * PSI(n - 1, m));

  auto const phi = [&coeffs](t_int n, t_int m) {
    return std::abs(m) > n ? 0 : coeffs(n * (n + 1) + m - 1);
  };
  auto const psi = [&coeffs, &N](t_int n, t_int m) {
    return std::abs(m) > n ? 0 : coeffs(n * (n + 1) + m - 1 + N * (N + 2));
  };
  t_complex const factor(0, 1 / (wavenumber * constant::c * constant::mu0));
  auto const factor_psi = factor * std::real(wavenumber * std::conj(wavenumber));
  auto const phiH = [&psi, factor_psi](t_int n, t_int m) { return factor_psi * psi(n, m); };
  auto const psiH = [&phi, factor](t_int n, t_int m) { return factor * phi(n, m); };

  for(t_int n(1), i(0); n <= N; ++n)
    for(t_int m(-n); m <= n; ++m, ++i) {
      E(i, 0) = OPTIMET_X(phi, psi);
      E(i, 1) = OPTIMET_Y(phi, psi);
      E(i, 2) = OPTIMET_Z(phi, psi);
      H(i, 0) = OPTIMET_X(phiH, psiH);
      H(i, 1) = OPTIMET_Y(phiH, psiH);
      H(i, 2) = OPTIMET_Z(phiH, psiH);
    }
}

template <class T>
Fields::Result Fields::single_sphere(Position const &position, t_int nMax,
                                     Eigen::MatrixBase<T> const &coeffs) const {
  Result result = Result::Zero();
  Eigen::Matrix<t_real, 3, 1> const spherical(
      position.stableNorm(), // r
      std::atan2(std::sqrt(position[0] * position[0] + position[1] * position[1]),
                 position[2]), // phi
      std::atan2(position[1], position[0]));
  for(t_int n(1), row(0); n <= nMax; ++n)
    for(t_int m(-n); m <= n; ++m, ++row) {
      auto const value =
          std::get<0>(optimet::bessel<Hankel1>(spherical(0) * wavenumber, n)).back() *
          boost::math::spherical_harmonic(n, m, spherical(1), spherical(2));
      result.col(0) += value * coeffs.row(row).head(3);
      result.col(1) += value * coeffs.row(row).tail(3);
    }
  return result;
}
}

#endif
