#ifndef BESSEL_H_
#define BESSEL_H_

#include <complex>
#include "constants.h"

using std::complex;

/**
 * The Bessel class implements the Spjerical Bessel and Hankel functions
 * and their derivatives, calculated from the zeroth order up to the maximum
 * order.
 * @warning The zeroth order derivative (never used) is not accurate!
 */
class Bessel {
private:
	complex<double> argument;	/**< The Bessel function argument. */
	long int besselType;				/**< The Bessel function type (0 - Bessel, 1 - Hankel (first kind), 2 - Hankel (second kind). */
	long int scale;					/**< Specifies if function should be scaled or not. */
	long int maxOrder;				/**< The maximum order to compute the function for. */

	bool initDone;				/**< Specify if the object has been initialized. */

public:

	complex<double> *data; 	/**< The direct Bessel functions. */
	complex<double> *ddata;	/**< The derivative Bessel functions. */

	long int zeroUnderflow;	/**< AMOS number of orders set to zero due to underflow. */
	long int ierr;			/**< AMOS or internal return error code. */

	/**
	 * Default constructor for the Bessel class.
	 * Does NOT initialize the object.
	 */
	Bessel();

	/**
	 * Initializing constructor for the Bessel class.
	 * @param argument_ the argument for the Bessel function.
	 * @param besselType_ the type of function (0 - Bessel, 1 - Hankel (first kind), 2 - Hankel (second kind).
	 * @param scale_  the scaling type (0 - unscaled, 1 - scaled).
	 * @param maxOrder_ the maximum order of functions to calculate.
	 */
	Bessel(complex<double> argument_, int besselType_, int scale_, int maxOrder_);

	/**
	 * Default destructor for the Bessel class.
	 */
	virtual ~Bessel();

	/**
	 * Initialization method for the Bessel class.
	 * @param argument_ the argument for the Bessel function.
	 * @param besselType_ the type of function (0 - Bessel, 1 - Hankel (first kind), 2 - Hankel (second kind).
	 * @param scale_  the scaling type (0 - unscaled, 1 - scaled).
	 * @param maxOrder_ the maximum order of functions to calculate.
	 * @return 0 if succesfull, value of ierr otherwise.
	 */
	int init(complex<double> argument_, int besselType_, int scale_, int maxOrder_);


	/**
	 * Calculate the Bessel functions and populate data and ddata.
	 * @return 0 if successful, value of ierr otherwise.
	 */
	int populate(void);
};

#endif /* BESSEL_H_ */
