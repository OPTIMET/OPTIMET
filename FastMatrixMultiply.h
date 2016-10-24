#ifndef FAST_MATRIX_MULTIPLY_H

#include "CoAxialTranslationCoefficients.h"
#include "RotationCoaxialDecomposition.h"
#include "RotationCoefficients.h"
#include "Scatterer.h"
#include "Types.h"
#include <utility>
#include <vector>

namespace optimet {
//! \brief Fast multiplication of effective incident field by transfer matrix
//! \details The multiplication occurs for a set of scattering particles receiving the EM radiation
//! to a set of locations where the fields are checked. Each incident effective field is
//! composed of coefficients for both Φ and Ψ potentials (in R basis).
class FastMatrixMultiply {
public:
  //! Range of incident/scattering objects
  typedef std::pair<t_int, t_int> Range;

  //! Creates the fast matrix multiply object
  //! \param[in] em_background: Electro-magnetic properties of the background material
  //! \param[in] wavenumber: Wave-number of the incident plane-wave
  //! \param[in] scatterers: all spherical scatterers in the problem
  //! \param[in] incident: Indices defining the scatterers this instance will take care of.
  //!     In practice, the input vector to the multiplication will contain the coefficients of the
  //!     effective incident field expanded on a spherical basis set at the location of the
  //!     scatterers in this range.
  //! \param[in] translate: Indices defining the location for which to expand the field on output.
  //!     In practice, the output vector of the multiplication will contain the coefficients of the
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

  //! Applies Mie coefficient to incident values
  Vector<t_complex> apply_mie_coefficients(Vector<t_complex> const &in) const;

  //! \brief Adds co-axial rotation/translation for a given particle pair
  //! \details the coefficients from input are projected from iScatt to oScatt using the given
  //! rotation translation functors. The result are added to the out-going coefficients.
  template <class T0, class T1>
  static void add_translation(Eigen::MatrixBase<T0> const &input, Scatterer const &iScatt,
                              Scatterer const &oScatt, Rotation const &rotation,
                              CachedCoAxialRecurrence::Functor const &translation,
                              t_real wavenumber, Eigen::MatrixBase<T1> &out);

  //! Adds co-axial rotation/translation for a given particle pair
  template <class T0, class T1, class T2>
  static void
  add_translation(Eigen::MatrixBase<T0> const &input, Scatterer const &iScatt,
                  Scatterer const &oScatt, Rotation const &rotation,
                  CachedCoAxialRecurrence::Functor const &translation, t_real wavenumber,
                  Eigen::MatrixBase<T1> &out, Eigen::MatrixBase<T2> &work);

  //! Apply translation to each particle pair
  void translation(Vector<t_complex> const &in, Vector<t_complex> &out) const;

protected:
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
  static constexpr t_int nfunctions(t_int nmax) { return (nmax + 1) * (nmax + 1); }
};

template <class T0, class T1>
void FastMatrixMultiply::add_translation(Eigen::MatrixBase<T0> const &input,
                                         Scatterer const &iScatt, Scatterer const &oScatt,
                                         Rotation const &rotation,
                                         CachedCoAxialRecurrence::Functor const &translation,
                                         t_real wavenumber, Eigen::MatrixBase<T1> &out) {
  if(input.cols() != 2)
    throw std::runtime_error("input should have two columns");
  if(out.cols() != 2)
    throw std::runtime_error("output should have two columns");

  auto const nmax = std::max(iScatt.nMax, oScatt.nMax);
  Eigen::Matrix<t_complex, Eigen::Dynamic, 4> work((nmax + 1) * (nmax + 1), 4);
  FastMatrixMultiply::add_translation(input, iScatt, oScatt, rotation, translation, out, work);
}

template <class T0, class T1, class T2>
void FastMatrixMultiply::add_translation(Eigen::MatrixBase<T0> const &input,
                                         Scatterer const &iScatt, Scatterer const &oScatt,
                                         Rotation const &rotation,
                                         CachedCoAxialRecurrence::Functor const &translation,
                                         t_real wavenumber, Eigen::MatrixBase<T1> &out,
                                         Eigen::MatrixBase<T2> &work) {
  auto const in_rows = nfunctions(iScatt.nMax);
  auto const out_rows = nfunctions(oScatt.nMax);

  auto const max_rows = std::max(in_rows, out_rows);
  if(work.rows() < max_rows or work.cols() != 4)
    work.resize(max_rows, 4);
  work.fill(0);

  assert((iScatt.vR.toEigenCartesian() - oScatt.vR.toEigenCartesian()).stableNorm() >=
         iScatt.radius + oScatt.radius);

  // First apply rotation
  auto rotated = work.block(0, 0, 2, in_rows);
  rotation.transpose(input, rotated);

  // Then perform co-axial translation
  auto ztrans = work.block(2, 0, 4, out_rows);
  translation(rotated, ztrans);

  // Then apply field-coaxial-tranlation transform thing
  auto const tz = (oScatt.vR.toEigenCartesian() - iScatt.vR.toEigenCartesian()).stableNorm();
  auto &field_translated = rotated;
  rotation_coaxial_decomposition(wavenumber, tz, ztrans, field_translated);

  // Finally, rotate back and adds to out-going coeffs
  auto &unrotated = ztrans;
  rotation.conjugate(field_translated, unrotated);
  out += unrotated.topRows(out_rows);
}
}

#endif
