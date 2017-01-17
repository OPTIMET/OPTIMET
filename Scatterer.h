#ifndef SCATTERER_H_
#define SCATTERER_H_

#include "ElectroMagnetic.h"
#include "Spherical.h"
#include "Types.h"
#include <vector>

/**
 * The Scatterer class is the highest level element of a geometry.
 * This is the only class used to create the geometry. The radius property
 * will define a virtual sphere encompassing the entire structure. Different
 * behaviors will be created based on the type of scatterer. Current
 * implementation includes only spherical type scatterers.
 */
class Scatterer {
public:
  /**
   * Default Scatterer constructor.
   * Does NOT initialize the object.
   */
  Scatterer(int nMax_);

  /**
   * Initialization constructor for Scatterer.
   * Creates a sphere scatterer centered in vR with properties elmag and
   * radius r.
   * @param vR_ the coordinates of the center of the scatterer.
   * @param elmag_ the electromagnetic properties of the scatterer.
   * @param radius_ the radius of the virtual sphere.
   * @param nMax_ the maximum value of the n iterator.
   * @see init()
   */
  Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_);

  /**
   * Initialization constructor for Scatterer.
   * Creates a sphere scatterer centered in vR with properties elmag and
   * radius r.
   * @param pos the position in cartesian coordinates
   * @param elmag_ the electromagnetic properties of the scatterer.
   * @param radius_ the radius of the virtual sphere.
   * @param nMax_ the maximum value of the n iterator.
   * @see init()
   */
  template <class T>
  Scatterer(Eigen::MatrixBase<T> const &pos, ElectroMagnetic elmag, double radius, int nMax)
      : Scatterer(Spherical<optimet::t_real>::toSpherical(
                      Cartesian<optimet::t_real>(pos[0], pos[1], pos[2])),
                  elmag, radius, nMax) {}

  /**
   * Default scatterer destructor.
   */
  virtual ~Scatterer();

  Spherical<double> vR;  /**< The coordinates of the center of the scatterer.*/
  ElectroMagnetic elmag; /**< The electromagnetic properties of the scatterer.*/
  double radius;         /**< The radius of a sphere encompassing the scatterer.*/
  int nMax;              /**< Maximum value of the n iterator. */

  std::vector<std::complex<double>> sourceCoef; /**< The source coefficients needed for SH work.*/

  /**
   * Returns the single object local scattering matrix T
   * Currently implements the scattering matrix for spherical objects
   * @param omega_ the angular frequency of the simulation
   * @param bground background electromagnetic medium
   * @return 0 if successful, 1 otherwise.
   */
  optimet::Vector<optimet::t_complex>
  getTLocal(optimet::t_real omega_, ElectroMagnetic const &bground) const;

  //! Coefficients for field inside a sphere
  optimet::Vector<optimet::t_complex>
  getIaux(optimet::t_real omega_, ElectroMagnetic const &bground) const;
};

#endif /* SCATTERER_H_ */
