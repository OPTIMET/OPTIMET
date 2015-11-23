#ifndef RESULT_H_
#define RESULT_H_

#include "SphericalP.h"
#include "Spherical.h"
#include "CompoundIterator.h"
#include "Geometry.h"
#include "OutputGrid.h"
#include "Excitation.h"

#include <complex>

using std::complex;

/**
 * The Result class is used to provide post simulation
 * output functions including field profiles, absorbption
 * and extinction cross sections, and scattering
 * coefficients.
 */

class Result
{
private:
	int nMax;				/**< Maximum number of harmonics. */
	Geometry *geometry;		/**< Pointer to the Geometry. */
	Excitation *excitation;	/**< Pointer to the Excitation. */
	complex<double> waveK;	/**< The complex wave number. */
	bool initDone;			/**< Specifies if object was initialized. */
	bool flagSH;			/**< Specifies if handling Second Harmonic results. */
	Result *result_FF;		/**< The Fundamental Frequency results vector. */
public:
	complex<double> *scatter_coef;		/**< The scattering coefficients. */
	complex<double> *internal_coef;		/**< The internal coefficients. */
	complex<double> *c_scatter_coef;	/**< The cluster centered scattering coefficients. */

	/**
	 * Default constructor for the Result class.
	 * Does NOT initialize the object.
	 */
	Result();

	/**
	 * Initialization constructor for the Result class.
	 * Fundamental Frequency version.
	 * @param geometry_ the pointer to the geometry.
	 * @param excitation_ the pointer to the excitation.
	 * @param nMax_ the maximum number of harmonics.
	 */
	Result(Geometry *geometry_, Excitation *excitation_, int nMax_);

	/**
	 * Initialization constructor for the Result class.
	 * Second Harmonic version.
	 * @param geometry_ the pointer to the geometry.
	 * @param excitation_ the pointer to the excitation.
	 * @param result_FF_ the pointer to the Fundamental Frequency results.
	 * @param nMax_ the maximum number of harmonics.
	 */
	Result(Geometry *geometry_, Excitation *excitation_, Result *result_FF_, int nMax_);

	/**
	 * Default destructor for the Result class.
	 */
	virtual ~Result();

	/**
	 * Initialization method for the Result class.
	 * Fundamental Frequency version.
	 * @param geometry_ the pointer to the geometry.
	 * @param excitation_ the pointer to the excitation.
	 * @param nMax_ the maximum number of harmonics.
	 */
	void init(Geometry *geometry_, Excitation *excitation_, int nMax_);

	/**
	 * Update method for the Result class.
	 * @param geometry_ the pointer to the geometry.
	 * @param excitation_ the pointer to the excitation.
	 * @param nMax_ the maximum number of harmonics.
	 */
	void update(Geometry *geometry_, Excitation *excitation_, int nMax_);

	/**
	 * Initialization constructor for the Result class.
	 * Second Harmonic version.
	 * @param geometry_ the pointer to the geometry.
	 * @param excitation_ the pointer to the excitation.
	 * @param result_FF_ the pointer to the Fundamental Frequency results.
	 * @param nMax_ the maximum number of harmonics.
	 */
	void init(Geometry *geometry_, Excitation *excitation_, Result *result_FF_, int nMax_);

	/**
	 * Returns the E field at a given point using the cluster centered formulation.
	 * @warning Testing method only. DO NOT USE IN PRODUCTION CODE!
	 * @param R_ the position of the point.
	 * @param projection_ defines spherical (1) or cartezian (0) projection.
	 * @return the value of the E field
	 */
	SphericalP<complex<double> > getEFieldC(Spherical<double> R_, int projection_);

	/**
	 * Returns the E and H fields at a given point.
	 * @param R_ the coordinates of the point.
	 * @param EField_ SphericalP vector that will store the E field.
	 * @param HField_ SphericalP vector that will store the H field.
	 * @param projection_ defines spherical (1) or cartezian (0) projection.
	 */
	void getEHFields(Spherical<double> R_, SphericalP<complex<double> > &EField_, SphericalP<complex<double> > &HField_, int projection_);

	/**
	 * Returns the E and H fields for a single harmonic and/or TE/TM component.
	 * @param R_ the coordinates of the point.
	 * @param EField_ SphericalP vector that will store the E field.
	 * @param HField_ SphericalP vector that will store the H field.
	 * @param projection_ defines spherical (1) or cartezian (0) projection.
	 * @param p_ the harmonic to be used.
	 * @param singleComponent_ return TE+TM (0), TE(1) or TM(2).
	 */
	void getEHFieldsModal(Spherical<double> R_, SphericalP<complex<double> > &EField_, SphericalP<complex<double> > &HField_, int projection_, CompoundIterator p_, int singleComponent_);

	/**
	 * Center the scattering coefficients.
	 * @warning Test method only! DO NOT USE FOR PRODUCTION CODE!
	 */
	void centerScattering();

	/**
	* Returns the Extinction Cross Section.
	* @return the extinction cross section.
	*/
	double getExtinctionCrossSection();

	/**
	* Returns the Absorption Cross Section.
	* @return the absorptions cross section.
	*/
	double getAbsorptionCrossSection();

	/**
	* Populate a grid with E and H fields.
	* @param oEGrid_ the OutputGrid object for the E fields.
	* @param oHGrid_ the OutputGrid object for the H fields.
	* @param projection_ defines spherical (1) or cartezian (0) projection.
	* @return 0 if succesful, 1 otherwise.
	*/
	int setFields(OutputGrid oEGrid_, OutputGrid oHGrid_, int projection_);

	/**
	* Populate a grid with E and H fields for a single harmonic and/or TE/TM component.
	* @param oEGrid_ the OutputGrid object for the E fields.
	* @param oHGrid_ the OutputGrid object for the H fields.
	* @param projection_ defines spherical (1) or cartezian (0) projection.
	* @param p_ the harmonic to be used.
	* @param singleComponent_ return TE+TM (0), TE(1) or TM(2).
	* @return 0 if succesful, 1 otherwise.
	*/
	int setFieldsModal(OutputGrid oEGrid_, OutputGrid oHGrid_, int projection_, CompoundIterator p_, int singleComponent_);

	/**
	 * Return the dominant harmonic.
	 * @return the CompoundIterator corresponding to the dominant harmonic (TE+TM).
	 */
	CompoundIterator getDominant();

	/**
	 * Returns the E and H fields at a given point using either the scattered or internal field methods.
	 * Used for continuity tests.
	 * @param R_ the coordinates of the point.
	 * @param EField_ SphericalP vector that will store the E field.
	 * @param HField_ SphericalP vector that will store the H field.
	 * @param projection_ defines spherical (1) or cartezian (0) projection.
	 * @param inside_ uses the internal (1) or external (0) field calculations.
	 */
	void getEHFieldsContCheck(Spherical<double> R_, SphericalP<complex<double> > &EField_, SphericalP<complex<double> > &HField_, int projection_, int inside_);

	/**
	 * Writes a set of files that check the field continuity around a particular object.
	 * @param objectIndex_ the object index for field continuity check.
	 */
	void writeContinuityCheck(int objectIndex_);
};

#endif /* RESULT_H_ */
