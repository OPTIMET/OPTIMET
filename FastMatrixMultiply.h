#ifndef FAST_MATRIX_MULTIPLY_H

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
            compute_mie_coefficients(em_background, wavenumber, scatterers_, incident_range_)) {}
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
  Eigen::Matrix<t_complex, Eigen::Dynamic, 2>
  apply_mie_coefficients(Vector<t_complex> const &in) const;

  //! Apply translation
  template <class T0, class T1>
  void apply_translation(Eigen::MatrixBase<T0> &out, Eigen::MatrixBase<T1> const &in) const {
    out = in;
  };

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
};
}

#endif
