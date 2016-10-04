#ifndef COAXIAL_TRANSLATION_COEFFICIENTS_H
#include "Types.h"
#include <array>
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

// class CoAxialTranslationAdditionCoefficients {
// public:
//   CoAxialTranslationAdditionCoefficients(t_real distance, t_complex waveK, bool regular = true)
//       : cached_recurrence(distance, waveK, regular) {}
//
//   //! \brief Computes the coefficients as per Gumerov (2003)
//   t_complex operator()(t_int n, t_int m, t_int l);
//
// protected:
//   //! Recurrence for all m
//   CachedCoAxialRecurrence cached_recurrence;
// };
}
#endif
