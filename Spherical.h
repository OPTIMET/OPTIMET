#ifndef SPHERICAL_H_
#define SPHERICAL_H_

#include "Cartesian.h"
#include "SphericalP.h"
#include <cmath>

using std::sqrt;
using std::sin;
using std::cos;
using std::acos;
using std::atan2;

/**
 * The Spherical class implements spherical coordinates.
 * This is a template class so any standard type may be used to create
 * spherical coordinates.
 * @warning Do not use without initialization.
 */
template<class sphType>
class Spherical
{
private:
  /**
   * Returns the point in Cartesian format.
   * @return the Cartesian coordinate object.
   */
  Cartesian<sphType> toCartesian()
  {
    return Cartesian<sphType>(rrr * sin(the) * cos(phi),
        rrr * sin(the) * sin(phi),
        rrr * cos(the));
  }

  /**
   * Returns the point in spherical projection format.
   * @return the SphericalP coordinate object.
   */
  SphericalP<sphType> toSphericalP()
  {
    return SphericalP<sphType>(rrr * sin(the) * cos(phi),
        rrr * sin(the) * sin(phi),
        rrr * cos(the));
  }

  /**
   * Converts a Cartesian vector into a Spherical vector.
   * @param point the Cartesian vector to be converted.
   * @return the vector in Spherical format.
   */
  Spherical<sphType> toSpherical(Cartesian <sphType> point)
  {
    double r_l = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    if(r_l > 0.0)
    {
      return Spherical<double>(r_l,
          acos(point.z / r_l),
          atan2(point.y, point.x));
    }

    return Spherical<double>(0.0, 0.0, 0.0);
  }


public:
  sphType rrr; /**< The R spherical coordinate. */
  sphType the; /**< The Theta spherical coordinate. */
  sphType phi; /**< The Phi spherical coordinates. */

  /**
   * Initialization constructor for Spherical.
   * @param rrr_ the sphType value for R
   * @param the_ the sphType value for Theta
   * @param phi_ the sphType value for Phi
   * @see init()
   */
  Spherical(sphType rrr_, sphType the_, sphType phi_)
  {
    init(rrr_, the_, phi_);
  }

  /**
   * Default Spherical constructor.
   */
  Spherical(void)
  {
    //
  }

  /**
   * Default Spherical destructor.
   */
  ~Spherical(void)
  {
    //
  }

  /**
   * Initialize a Spherical object.
   * @param rrr_ the sphType value for R
   * @param the_ the sphType value for Theta
   * @param phi_ the sphType value for Phi
   * @see Spherical()
   */
  void init(sphType rrr_, sphType the_, sphType phi_)
  {
    rrr = rrr_;
    the = the_;
    phi = phi_;
  }

  /**
   * Implements the dot-product of two spherical vectors.
   * @param argument_ the vector to be multiplied with.
   * @return the dot-product of this and argument_.
   */
  sphType operator * (Spherical<sphType> argument_)
  {
    return toCartesian() * argument_.toCartesian();
  }

  /**
   * Implements the dot-product of a spherical and cartesian vector.
   * @param argument_ the cartesian vector to be multiplied with.
   * @return
   */
  sphType operator * (Cartesian<sphType> argument_)
  {
    return toCartesian() * argument_;
  }

  /**
   * Implements the dot-product of a spherical and spherical projection vector.
   * @param argument_ the SphericalP vector to be multiplied with.
   * @return
   */
  sphType operator * (SphericalP<sphType> argument_)
  {
    return toSphericalP() * argument_;
  }

  /**
   * Implements the vector*scalar product.
   * @param argument_ the scalar to be multiplied by.
   * @return the vector*scalar product.
   */
  Spherical<sphType> operator * (sphType argument_)
  {
    return Spherical<sphType>(rrr * argument_, the * argument_, phi * argument_);
  }

  /**
   * Implements the difference between this vector and the argument_.
   * @param argument_ the vector to be subtracted.
   * @return a Spherical vector as the difference.
   */
  Spherical<sphType> operator - (Spherical<sphType> argument_)
  {
    return toSpherical(toCartesian() - argument_.toCartesian());
  }
};


#endif /* SPHERICAL_H_ */
