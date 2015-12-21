#ifndef TRANSLATION_ADDITION_COEFFICIENTS_H

#include "Types.h"
#include <map>
#include <array>

#include "Spherical.h"

namespace optimet {
//! Equation 1 of Appendix A in Stout (2002)
t_complex Ynm(Spherical<t_real> const &R, t_int n, t_int m);

//! \brief Computes translation addition coefficients
//! \details Equations come from Stout (2002), appendix C.
class TranslationAdditionCoefficients {
  //! Indices tuple
  typedef std::array<t_int, 4> t_indices;

public:
  TranslationAdditionCoefficients(Spherical<t_real> R, t_complex waveK, bool regular = true)
      : direction(R), waveK(waveK), regular(regular) {}

  //! \brief Returns translation addition coefficients
  //! \details n and m correspond to the same variables in Stout (2004), l and k correspond to ν and
  //! μ, respectively.
  t_complex operator()(t_int n, t_int m, t_int l, t_int k);
  //! Returns coefficient
  t_complex operator()(t_indices const &indices) {
    return (*this)(indices[0], indices[1], indices[2], indices[3]);
  }

protected:
  //! Direction of the incident wave
  Spherical<t_real> direction;
  //! Wavenumber of the incident wave
  t_complex waveK;
  //! Whether this is for regular or irregular coeffs
  bool regular;
  //! Caches known coefficients
  std::map<t_indices, t_complex> cache;
  //! Caches known coefficients for m < 0 version
  std::map<t_indices, t_complex> cache_conj;

  //! Initial values
  t_complex initial(t_int l, t_int k);
};
}

#endif
