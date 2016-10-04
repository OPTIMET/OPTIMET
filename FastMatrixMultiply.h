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

  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Range incident, Range translate)
      : em_background_(em_background), wavenumber_(wavenumber), scatterers_(scatterers),
        global_indices_(compute_indices(scatterers_)), incident_range_(incident),
        translate_range_(translate),
        rotations_(compute_rotations(scatterers_, incident_range_, translate_range_)),
        mie_coefficients_(
            compute_mie_coefficients(em_background, wavenumber, scatterers_, incident_range_)) {}
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers)
      : FastMatrixMultiply(em_background, wavenumber, scatterers, {0, scatterers.size()},
                           {0, scatterers.size()}) {}

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

  //! Applies Mie coefficient to incident values
  Vector<t_complex> apply_mie_coefficients(Vector<t_complex> const &in) const;

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
