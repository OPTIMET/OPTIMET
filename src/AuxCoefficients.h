// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#ifndef AUX_COEFFICIENTS_H_
#define AUX_COEFFICIENTS_H_

#include "Spherical.h"
#include "SphericalP.h"
#include "Bessel.h"
#include "Types.h"

#include <complex>
#include <vector>
#include <tuple>

namespace optimet {

/**
 * The AuxCoefficients class implements the spherical functions M, N, B and C.
 * Also implements the dn symbol for incoming wave definition.
 */
class AuxCoefficients {
public:
  /**
   * Compute the d_n symbol.
   * @param nMax the maximum value of n iterator.
   * @return the d_n symbol vector.
   */
  static std::vector<t_real> compute_dn(t_uint nMax);

  /**
   * Compute the Wigner functions and their derivatives.
   * @param nMax the maximum value of the n iterator.
   * @param m the value of the m iterator.
   * @param R the Spherical vector.
   * @return the Wigner and dWigner vectors, in that order.
   */
  static std::tuple<std::vector<t_real>, std::vector<t_real>>
  VIGdVIG(t_uint nMax, t_int m, const Spherical<t_real> &R);

  /**
   * Initializing constructor for the AuxCoefficients class.
   * @param R_ the Spherical vector.
   * @param waveK_ the wave number.
   * @param regular_ the type of coefficients (regular or not).
   * @param nMax_ the maximum value of the n iterator.
   */
  AuxCoefficients(const Spherical<t_real> &R, t_complex waveK, bool regular,
                  t_uint nMax);

  const SphericalP<t_complex> &M(t_uint i) const { return _M[i]; }
  const SphericalP<t_complex> &Xp(t_uint i) const { return _Xp[i]; }
  const SphericalP<t_complex> &Xm(t_uint i) const { return _Xm[i]; }
  const SphericalP<t_complex> &N(t_uint i) const { return _N[i]; }
  const SphericalP<t_complex> &B(t_uint i) const { return _B[i]; }
  const SphericalP<t_complex> &C(t_uint i) const { return _C[i]; }
  
  
  const t_real &dn(t_uint i) const { return _dn[i]; }

private:
  /**
   * Compute the Pn functions
   * @param nMax the maximum value of the n iterator.
   * @param Wigner the Wigner functions.
   * @return Pn.
   */
  std::vector<SphericalP<t_complex>>
  compute_Pn(t_uint nMax, const std::vector<t_real> &Wigner);

  /**
   * Compute the Cn functions
   * @param nMax the maximum value of the n iterator.
   * @param m_ the the value of the m iterator.
   * @param R the Spherical vector.
   * @param dn the d_n symbols.
   * @param Wigner the Wigner functions.
   * @param dWigner the derivatives of the Wigner functions.
   * @param Cn the vector to store Cn in.
   * @return 0 if succesful, 1 otherwise (deprecated).
   */
  std::vector<SphericalP<t_complex>>
  compute_Cn(t_uint nMax, t_int m, const Spherical<t_real> &R,
             const std::vector<t_real> &Wigner,
             const std::vector<t_real> &dWigner);
                     

  /**
   * Compute the Bn functions
   * @param nMax the maximum value of the n iterator.
   * @param m_ the the value of the m iterator.
   * @param R the Spherical vector.
   * @param dn the d_n symbols.
   * @param Wigner the Wigner functions.
   * @param dWigner the derivatives of the Wigner functions.
   * @param Bn the vector to store Bn in.
   * @return 0 if succesful, 1 otherwise (deprecated).
   */
  std::vector<SphericalP<t_complex>>
  compute_Bn(t_uint nMax, t_int m, const Spherical<t_real> &R,
             const std::vector<t_real> &Wigner,
             const std::vector<t_real> &dWigner);

  /**
   * Compute the Mn functions.
   * @param nMax the maximum value of the n iterator.
   * @param R the Spherical vector.
   * @param waveK the wave vector.
   * @param Cnm the Cnm functions.
   * @param Mn the Mn functions.
   * @param BHreg the type of Bessel functions (regular or not).
   * @return 0 if succesful, 1 otherwise (deprecated).
   */
  std::vector<SphericalP<t_complex>>
  compute_Mn(t_uint nMax, t_int m, const Spherical<t_real> &R,
             std::complex<t_real> waveK, const std::vector<t_real> &dn,
             const std::vector<SphericalP<t_complex>> &Cnm, BESSEL_TYPE BHreg);
             
        // compute the Xn functions       
             
    std::vector<SphericalP<t_complex>>
   compute_Xm1 (t_uint nMax, t_int m, const Spherical<t_real> &R,
             std::complex<t_real> waveK, const std::vector<t_real> &dn,
             const std::vector<SphericalP<t_complex>> &Pn, BESSEL_TYPE BHreg);
             
   std::vector<SphericalP<t_complex>>
   compute_Xp1 (t_uint nMax, t_int m, const Spherical<t_real> &R,
             std::complex<t_real> waveK, const std::vector<t_real> &dn,
             const std::vector<SphericalP<t_complex>> &Bn, BESSEL_TYPE BHreg); 

  /**
   * Compute the Nn functions.
   * @param nMax the maximum value of the n iterator.
   * @param R the Spherical vector.
   * @param waveK the wave vector.
   * @param Pn the Pn functions.
   * @param Bn the Bn functions.
   * @param Nn the Nn functions.
   * @param BHreg the type of Bessel functions (regular or not).
   * @return 0 if successful, 1 otherwise(deprecated).
   */
  std::vector<SphericalP<t_complex>>
  compute_Nn(t_uint nMax, t_int m, const Spherical<t_real> &R, t_complex waveK,
             const std::vector<t_real> &dn,
             const std::vector<SphericalP<t_complex>> &Pn,
             const std::vector<SphericalP<t_complex>> &Bn, BESSEL_TYPE BHreg);

  std::vector<SphericalP<t_complex>> _M, _N, _B,
      _C, _Xp, _Xm; /**< The M, N, B and C functions in compound iterator format. */

  std::vector<t_real>
      _dn; /**< The dn symbols (required for the Excitation class). */
};

} // namespace optimet

#endif /*AUX_COEFFICIENTS_H_*/
