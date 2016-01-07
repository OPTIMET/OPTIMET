#ifndef AUX_COEFFICIENTS_H_
#define AUX_COEFFICIENTS_H_

#include "Spherical.h"
#include "SphericalP.h"

#include <complex>

/**
 * The AuxCoefficients class implements the spherical functions M, N, B and C.
 * Also implements the dn symbol for incoming wave definition.
 * @author Ahmed Al-Jarro (original code)
 * @author Claudiu Biris (OOP framework implementation)
 */
class AuxCoefficients {
private:
  //  /**
  //   * Compute the d_n symbol.
  //   * @param nMax the maximum value of n iterator.
  //   * @param dn the vector to store dn in.
  //   * @return 0 if succesful, 1 otherwise (deprecated).
  //   */
  //  int compute_dn(int nMax, double *dn);

  /**
   * Compute the Pn functions
   * @param nMax the maximum value of the n iterator.
   * @param m_ the the value of the m iterator.
   * @param R the Spherical vector.
   * @param dn the d_n symbols.
   * @param Wigner the Wigner functions.
   * @param dWigner the derivatives of the Wigner functions.
   * @param Pn the vector to store Pn in.
   * @return 0 if succesful, 1 otherwise (deprecated).
   */
  int compute_Pn(int nMax, double *Wigner,
                 SphericalP<std::complex<double>> *Pn);

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
  int compute_Cn(int nMax, int m_, Spherical<double> R, double *Wigner,
                 double *dWigner, SphericalP<std::complex<double>> *Cn);

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
  int compute_Bn(int nMax, int m_, Spherical<double> R, double *Wigner,
                 double *dWigner, SphericalP<std::complex<double>> *Bn);

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
  int compute_Mn(int nMax, int m_, Spherical<double> R,
                 std::complex<double> waveK, double *dn,
                 SphericalP<std::complex<double>> *Cnm,
                 SphericalP<std::complex<double>> *Mn, int BHreg);

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
  int compute_Nn(int nMax, int m_, Spherical<double> R,
                 std::complex<double> waveK, double *dn,
                 SphericalP<std::complex<double>> *Pn,
                 SphericalP<std::complex<double>> *Bn,
                 SphericalP<std::complex<double>> *Nn, int BHreg);

  /**
   * Compute the M and N functions in compound iterator format.
   * @param R the Spherical vector.
   * @param waveK the wave number.
   * @param BHreg the type of Bessel function (regular or not).
   * @param nMax the maximum value of the n iterator.
   * @param dataMp the dataMp vector.
   * @param dataNp the dataNp vector.
   * @return
   */
  int compute_MpNp(Spherical<double> R, std::complex<double> waveK, int BHreg,
                   int nMax, SphericalP<std::complex<double>> *dataMp,
                   SphericalP<std::complex<double>> *dataNp);

  //  /**
  //   * Compute the Wigner functions and their derivatives.
  //   * @param nMax the maximum value of the n iterator.
  //   * @param m_ the value of the m iterator.
  //   * @param R the Spherical vector.
  //   * @param Wigner the Wigner vector.
  //   * @param dWigner the derivative Wigner vector.
  //   * @return 0 if succesful, 1 otherwise.
  //   */
  //  int VIGdVIG(int nMax, int m_, Spherical<double> R, double *Wigner, double
  //  *dWigner);

public:
  SphericalP<std::complex<double>>
      *dataMp; /**< The M functions in compound iterator format. */
  SphericalP<std::complex<double>>
      *dataNp; /**< The N functions in compound iterator format. */
  SphericalP<std::complex<double>>
      *dataBp; /**< The B functions in compound iterator format. */
  SphericalP<std::complex<double>>
      *dataCp; /**< The C functions in compound iterator format. */

  double *dn; /**< The dn symbols (required for the Excitation class). */

  Spherical<double> R; /**< The Spherical vector. */
  int besselType;      /**< Specifies regular or non-regular functions. */
  std::complex<double> waveK; /**< The complex wave number. */
  int nMax;                   /**< The maximum value of the n iterator. */

  /**
   * Initializing constructor for the AuxCoefficients class.l
   * @param R_ the Spherical vector.
   * @param waveK_ the wave number.
   * @param regular_ the type of coefficients (regular or not).
   * @param nMax_ the maximum value of the n iterator.
   */
  AuxCoefficients(Spherical<double> R_, std::complex<double> waveK_,
                  int regular_, int nMax_);

  /**
   * Compute the P_nm functions in compound iterator format.
   * Needed for transfer coefficients in Coupling object.
   * @param R the r coordinate.
   * @param nMax the maximum value of the n iterator.
   * @param dataPp the Pp vector.
   * @return 0 if successful, 1 otherwise.
   */
  int compute_Pp(Spherical<double> R, int nMax,
                 SphericalP<std::complex<double>> *dataPp);

  // AJ
  // -----------------------------------------------------------------------------
  // temporary call for checking addition-translation coefficients evaluations
  static int compute_dn(int nMax, double *dn);
  static int VIGdVIG(int nMax, int m_, Spherical<double> R, double *Wigner,
                     double *dWigner);
};

#endif /*AUX_COEFFICIENTS_H_*/
