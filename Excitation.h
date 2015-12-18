#ifndef EXCITATION_H_
#define EXCITATION_H_

#include "Spherical.h"

#include <complex>

/**
 * The Excitation class implements the incoming wave coefficients.
 * Possible Excitation types are:
 *  0 - plane wave
 * The wave is considered to be in freespace. Do not forget to divide by
 * relative
 * speed of light in background medium.
 * @warning Do not use without initialization.
 */
class Excitation {
private:
  bool initDone; /**< Specifies if the object has been initialized. */
public:
  //! The incoming wave values as (0, Einc_the, Einc_phi) in [V/m].
  SphericalP<std::complex<double>> Einc;
  Spherical<double> vKInc; /**< The incoming wavevector angular values as (R,
                              Kinc_the, Kinc_pphi). */
  long nMax;               /**< The maximum value for the n iterator. */
  unsigned long type;      /**< The Excitation type. */

  std::complex<double> *dataIncAp; /**< The incoming wave a_n^m coefficients. */
  std::complex<double> *dataIncBp; /**< The incoming wave b_n^m coefficients. */

  std::complex<double> waveK; /**< The incoming wave wavenumber. */
  double lambda;              /**< The incoming wave wavelength (calculated). */
  double omega;               /**< The incoming wave frequency (calculated). */

  /**
   * Default constructor for the Excitation class.
   * Does NOT initialize the object.
   */
  Excitation();

  /**
   * Initialization constructor for the Excitation class.
   * @param type_ the type of the wave.
   * @param Einc_ the incoming wave values.
   * @param waveKInc_ the incoming wavevector values.
   * @param nMax_ the maximum value of the n iterator.
   */
  Excitation(unsigned long type_, SphericalP<std::complex<double>> Einc_,
             Spherical<double> vKInc_, int nMax_);

  /**
   * Default destructor for the Excitation class.
   */
  virtual ~Excitation();

  /**
   * Initialization method for the Excitation class.
   * @param type_ the type of the wave.
   * @param Einc_ the incoming wave values.
   * @param waveKInc_ the incoming wavevector values.
   * @param nMax_ the maximum value of the n iterator.
   */
  void init(unsigned long type_, SphericalP<std::complex<double>> Einc_,
            Spherical<double> vKInc_, int nMax_);

  /**
   * Update method for the Excitation class.
   * @param type_ the type of the wave.
   * @param Einc_ the incoming wave values.
   * @param waveKInc_ the incoming wavevector values.
   * @param nMax_ the maximum value of the n iterator.
   */
  void update(unsigned long type_, SphericalP<std::complex<double>> Einc_,
              Spherical<double> vKInc_, int nMax_);

  /**
   * Populate the a_n^m and b_n^m coefficients for the incoming wave.
   * Stored in dataIncAp and dataIncBp with CompoundIterator indexing.
   * @return 0 if successful, 1 otherwise.
   */
  int populate();

  /**
   * Returns the Incoming Local matrix converting the a,b coefficients into
   * local ones.
   * @param point_ the new origin of the coordinate system.
   * @param Inc_local_ the incoming local matrix.
   * @param nMax_ the maximum value of the n iterator.
   * @return 0 if successful, 1 otherwise.
   */
  int getIncLocal(Spherical<double> point_, std::complex<double> *Inc_local_,
                  int nMax_);

  /**
   * Updates the wavelength of the current excitation object to a new value.
   * @param lambda_ the new value of the wavelength.
   */
  void updateWavelength(double lambda_);
};

#endif /* EXCITATION_H_ */
