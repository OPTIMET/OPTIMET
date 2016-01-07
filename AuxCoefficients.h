#ifndef AUX_COEFFICIENTS_H_
#define AUX_COEFFICIENTS_H_

#include "Spherical.h"
#include "SphericalP.h"
#include "Bessel.h"

#include <complex>
#include <vector>
#include <tuple>

namespace optimet {

/**
 * The AuxCoefficients class implements the spherical functions M, N, B and C.
 * Also implements the dn symbol for incoming wave definition.
 * @author Ahmed Al-Jarro (original code)
 * @author Claudiu Biris (OOP framework implementation)
 */
class AuxCoefficients {
public:
  /**
   * Compute the d_n symbol.
   * @param nMax the maximum value of n iterator.
   * @return the d_n symbol vector.
   */
  static std::vector<double> compute_dn(int nMax);

  /**
   * Compute the Wigner functions and their derivatives.
   * @param nMax the maximum value of the n iterator.
   * @param m the value of the m iterator.
   * @param R the Spherical vector.
   * @return the Wigner and dWigner vectors, in that order.
   */
  static std::tuple<std::vector<double>, std::vector<double>>
  VIGdVIG(int nMax, int m, const Spherical<double> &R);

  /**
   * Initializing constructor for the AuxCoefficients class.
   * @param R_ the Spherical vector.
   * @param waveK_ the wave number.
   * @param regular_ the type of coefficients (regular or not).
   * @param nMax_ the maximum value of the n iterator.
   */
  AuxCoefficients(const Spherical<double> &R, std::complex<double> waveK,
                  bool regular, int nMax);

  const SphericalP<std::complex<double>> &M(std::size_t i) const {
    return _M[i];
  }
  const SphericalP<std::complex<double>> &N(std::size_t i) const {
    return _N[i];
  }
  const SphericalP<std::complex<double>> &B(std::size_t i) const {
    return _B[i];
  }
  const SphericalP<std::complex<double>> &C(std::size_t i) const {
    return _C[i];
  }

  const double &dn(std::size_t i) const { return _dn[i]; }

private:
  /**
   * Compute the Pn functions
   * @param nMax the maximum value of the n iterator.
   * @param Wigner the Wigner functions.
   * @return Pn.
   */
  std::vector<SphericalP<std::complex<double>>>
  compute_Pn(int nMax, const std::vector<double> &Wigner);

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
  std::vector<SphericalP<std::complex<double>>>
  compute_Cn(int nMax, int m, const Spherical<double> &R,
             const std::vector<double> &Wigner,
             const std::vector<double> &dWigner);

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
  std::vector<SphericalP<std::complex<double>>>
  compute_Bn(int nMax, int m, const Spherical<double> &R,
             const std::vector<double> &Wigner,
             const std::vector<double> &dWigner);

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
  std::vector<SphericalP<std::complex<double>>>
  compute_Mn(int nMax, int m, const Spherical<double> &R,
             std::complex<double> waveK, const std::vector<double> &dn,
             const std::vector<SphericalP<std::complex<double>>> &Cnm,
             BESSEL_TYPE BHreg);

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
  std::vector<SphericalP<std::complex<double>>>
  compute_Nn(int nMax, int m, const Spherical<double> &R,
             std::complex<double> waveK, const std::vector<double> &dn,
             const std::vector<SphericalP<std::complex<double>>> &Pn,
             const std::vector<SphericalP<std::complex<double>>> &Bn,
             BESSEL_TYPE BHreg);

  std::vector<SphericalP<std::complex<double>>> _M, _N, _B,
      _C; /**< The M, N, B and C functions in compound iterator format. */

  std::vector<double>
      _dn; /**< The dn symbols (required for the Excitation class). */
};

} // namespace optimet

#endif /*AUX_COEFFICIENTS_H_*/
