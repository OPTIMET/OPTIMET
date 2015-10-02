#ifndef LEGENDRE_H_
#define LEGENDRE_H_

#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gamma.h"
#include "CompoundIterator.h"
#include "Tools.h"
#include <iostream>
#include <cmath>


/**
 * The Legendre class implements associated Legendre polynomials using gsl.
 * The Legendre object is initialized and returns an array of Legendre
 * polynomials using a CompoundIterator. The result is of the type
 * \f$P_p(argument) = P_n^{-m} (argument)\f$ where the mapping
 * \f$p \rightarrow (n, m)\f$ is used.
 * @warning Do not use without initialization.
 * @warning Due to factorial sizes, Legendre has been tested up to maxOrder=6559
 * corresponding to p -> (80, -80), however, some of the return types are
 * inf. Overflow is triggered by GSL at around maxOrder=7500.
 * Legendre is released under the GSL. To view a copy
 * of the licence, look in the documentation.
 */

class Legendre
{
private:
	bool initDone;				/**< Verifies initialization. */
public:
	double *data;				/**< Associated Legendre polynomials. */
	double *ddata;				/**< Derivatives of Associated Legendre polynomials. */

	double argument;			/**< Legendre argument. */
	int nMax;					/**< Maximum order of n iterator. */

	/**
	 * Default constructor for the Legendre class.
	 * Does NOT intialize the object.
	 */
	Legendre();

	/**
	 * Initialization constructor for the Legendre class.
	 * @param argument_ the Legendre argument.
	 * @param maxOrder_ the maximum order of the compound iterator.
	 */
	Legendre(double argument_, int nMax_);

	/**
	 * Default destructor for the Legendre class.
	 */
	virtual ~Legendre();

	/**
	 * Initialization method for the Legendre class.
	 * @param argument_ the Legendre argument.
	 * @param maxOrder_ the maximum order of the compound iterator.
	 */
	void init(double argument_, int nMax_);

	/**
	 * Calculates the associated Legendre polynomials and populates data.
	 * The Legendre polynomials are calculated as
	 * \f$P_p(argument) = P_n^{-m} (argument)\f$ where the mapping
	 *\f$p \rightarrow (n, m)\f$ is used.
	 * @return 0 if succesful, 1 otherwise.
	 */
	int populate();
};

#endif /* LEGENDRE_H_ */
