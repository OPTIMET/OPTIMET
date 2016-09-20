#ifndef COAXIAL_TRANSLATION_COEFFICIENTS_H
#include "Types.h"
#include <array>
#include <map>

#include "Spherical.h"

namespace optimet {

namespace details {

class CachedCoAxialRecurrence {
public:
  //! Indices tuple
  typedef std::array<t_int, 3> t_indices;

  CachedCoAxialRecurrence(t_real distance, t_complex waveK, bool regular = true)
      : distance(distance), waveK(waveK), regular(regular) {}

  //! \brief Returns coaxial translation coefficients
  //! \details n, l and m correspond to the same variables in Gumerov (2002),
  //! s = m by definition.
  t_complex operator()(t_int n, t_int m, t_int l);

  bool is_regular() const { return regular; }

protected:
  //! Distance that the solution is to be translated by
  t_real const distance;
  //! Wavenumber of the incident wave
  t_complex const waveK;
  //! Whether this is for regular or irregular coeffs
  bool const regular;
  //! Caches known coefficients
  std::map<t_indices, t_complex> cache;

  //! Switches between recurrence relationships
  t_complex recurrence(t_int n, t_int m, t_int l);
  t_complex initial(t_int l);
  t_complex sectorial_recurrence(t_int n, t_int m, t_int l);
  t_complex zonal_recurrence(t_int n, t_int l);
  t_complex offdiagonal_recurrence(t_int n, t_int m, t_int l);
};
}
class CoAxialTranslationAdditionCoefficients {
public:
  CoAxialTranslationAdditionCoefficients(t_real distance, t_complex waveK, bool regular = true)
      : cached_recurrence(distance, waveK, regular) {}

  //! \brief Computes the coefficients as per Gumerov (2003)
  t_complex operator()(t_int n, t_int m, t_int l);

protected:
  //! Recurrence for all m
  details::CachedCoAxialRecurrence cached_recurrence;
};
}
#endif
