#ifndef COUPLING_H_
#define COUPLING_H_

#include "Spherical.h"

#include <complex>

/**
 * The Coupling class implements the calculation of the A and B coupling
 * coefficients.
 * @warning Do not use without initialization!
 */
class Coupling {
private:
  /**
   * Return the a_nm_p symbol.
   * @param n the value of n.
   * @param m the value of m.
   * @return a_nm_p.
   */
  double a_nm_p(double n, double m);

  /**
   * Return the a_nm_m symbol.
   * @param n the value of n.
   * @param m the value of m.
   * @return a_nm_m.
   */
  double a_nm_m(double n, double m);

  /**
   * Return the b_nm_p symbol.
   * @param n the value of n.
   * @param m the value of m.
   * @return b_nm_p.
   */
  double b_nm_p(double n, double m);

  /**
   * Return the b_nm_m symbol.
   * @param n the value of n.
   * @param m the value of m.
   * @return b_nm_m.
   */
  double b_nm_m(double n, double m);

  /**
   * Calculates the coupling coefficients for a relative vector R.
   * @param R the relative Spherical vector.
   * @param waveK the complex wave vector.
   * @param BHreg the type of Bessel function.
   * @param n_max the maximum value of the n iterator.
   * @param dataApq the storage vector for dataApq.
   * @param dataBpq the storage vector for dataBpq.
   */
  void TransferCoefficients(Spherical<double> R, std::complex<double> waveK,
                            int BHreg, int n_max,
                            std::complex<double> **dataApq,
                            std::complex<double> **dataBpq);

  bool initDone; /**< Specifies if the object has been initialized. */
public:
  long nMax;              /**< The maximum value of the n iterator. */
  Spherical<double> relR; /**< Relative spherical vector between two spheres. */
  std::complex<double> waveK; /**< The complex wave number. */

  std::complex<double> **dataApq; /**< The A_nmkl coefficients. */
  std::complex<double> **dataBpq; /**< The B_nmkl coefficients. */

  int regular; /**< Specifies if we should be using the regular (J) functions.
                  1-yes, 0-no. */

  /**
   * Default constructor for the Coupling class.
   * Does NOT initialize the object.
   */
  Coupling();

  /**
   * Initialization constructor for the Coupling class.
   * @param relR_ the relative spherical vector between two spheres.
   * @param waveK_ the complex wave number.
   * @param regular_ the regular flag.
   * @param nMax_ the maximum value of the n iterator.
   */
  Coupling(Spherical<double> relR_, std::complex<double> waveK_, int regular_,
           long nMax_);

  /**
   * Default destructor for the Coupling class.
   */
  virtual ~Coupling();

  /**
   * Initialization method for the Coupling class.
   * @param relR_ the relative spherical vector between two spheres.
   * @param waveK_ the complex wave number.
   * @param regular_ the regular flag.
   * @param nMax_ the maximum value of the n iterator.
   */
  void init(Spherical<double> relR_, std::complex<double> waveK_, int regular_,
            long nMax_);

  /**
   * Populate the dataApq and dataBpq variables.
   * @return 0 if successful, 1 otherwise.
   */
  int populate();
};

#endif /* COUPLING_H_ */
