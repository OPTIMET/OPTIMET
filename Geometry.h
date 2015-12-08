#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include "Scatterer.h"
#include "Excitation.h"

/**
 * The Geometry class implements a list of objects and properties of the medium.
 * To use the class, first initialize with init or the constructor by
 * declaring an expected capacity. Add objects with pushObject (Scatterer type)
 * until geometry is full. Finally, always validate to make sure the number
 * of objects are properly stored and objects do not intersect.
 * @warning Never use this class without initializing and validating.
 */
class Geometry
{
private:
  //
public:
  std::vector<Scatterer> objects;     /**< The list of scatterers. */

  ElectroMagnetic bground;  /**< The properties of the background. */

  bool validDone;       /**< Specifies if the geometry has been validated. */

  int structureType;      /**< The type of structure used: 0 - non-structured, 1 - spiral . */

  double spiralSeparation;  /**< Spiral separation (NOT distance between centers but between surfaces. */
  int normalToSpiral;     /**< Spiral normal plane (0 - x, 1 - y, 2 - z). */

  /**
   * Default constructor for the Geometry class. Does not initialize.
   */
  Geometry();

  /**
   * Alternative constructor for the Geometry class.
   * Can initialize with a list of Scatterer type objects.
   * @param objects_ the list of Scatterer objects
   * @see init()
   * @warning Try to use the pushObject method for best results.
   * @see pushObject()
   */
  // Geometry(int capacity_, Scatterer* objects_);

  /**
   * Default destructor for the Geometry class.
   */
  virtual ~Geometry();

  /**
   * Initialize the background. Default is vacuum.
   * @param bground_ the ElectroMagnetic properties of the background.
   */
  void initBground(ElectroMagnetic bground_);

  /**
   * Add an object to the Geometry.
   * Will discard any objects pushed beyond the declared capacity.
   * @param object_ the Scatterer type object to be added.
   * @return 0 if add successful, 1 otherwise
   * @see init()
   * @see noObjects
   * @see capacity
   */
  int pushObject(Scatterer object_);

  /**
   * Validate the geometry.
   * @return 1 if geometry is valid, 0 otherwise.
   */
  int validate();

  /**
   * Returns the single object local scattering matrix T_j.
   * Currently implements the scattering matrix for spherical objects.
   * @param omega_ the angular frequency of the simulation.
   * @param objectIndex_ the index of the object for which T_j is calculated.
   * @param nMax_ the maximum value of the n iterator.
   * @param T_local_ the return value as the local T_j scattering matrix.
   * @return 0 if successful, 1 otherwise.
   */
  int getTLocal(double omega_, int objectIndex_, int nMax_, std::complex<double> ** T_local_);

  int getIaux(double omega_, int objectIndex_, int nMax_, std::complex<double> *I_aux_);

  int getCabsAux (double omega_, int objectIndex_, int nMax_, double *Cabs_aux_);

  int getNLSources(double omega_, int objectIndex_, int nMax_, std::complex<double> *sourceU, std::complex<double> *sourceV);

  /**
   * Returns the relative vector R_lj between two objects.
   * @param firstObject_ the index of the first (l) object.
   * @param secondObjecT_ the index of the second (j) object.
   * @return the vector R_lj = R_l-R_j.
   */
  Spherical<double> translateCoordinates(int firstObject_, int secondObject_);

  /**
   * Checks if a point R is located inside an object.
   * @param R_ coordinates of the point to check.
   * @return the index of the object in which this point is or -1 if outside.
   */
  int checkInner(Spherical<double> R_);

  /**
   * Calculates the local second harmonic sources.
   * @param objectIndex_ the index of the object.
   * @param incWave_ pointer to the incoming excitation.
   * @param scatterCoef_ pointer to the ENTIRE scattering coefficients for the FF case
   * @param nMax_ the maximum value of the n iterator.
   * @param Q_SH_local_ the return value of the local SH source vector.
   * @return 0 if successful, 1 otherwise.
   */
  int getSourceLocal(int objectIndex_, Excitation *incWave_, std::complex<double> *internalCoef_FF_, int nMax_, std::complex<double>* Q_SH_local_);

  int setSourcesSingle(Excitation *incWave_, std::complex<double> *internalCoef_FF_, int nMax_);

  /**
   * Updates the Geometry object to a new Excitation.
   * @param lambda_ the new wavelength.
   */
  void update(Excitation *incWave_);

  /**
   * Updates the Geometry object by modifying the radius of an object.
   * Also validates the geometry.
   * @param radius_ the new radius.
   * @param object_ the object to be modified.
   */
  void updateRadius(double radius_, int object_);

  /**
   * Rebuilds a structure based on the new updated Radius.
   */
  void rebuildStructure();
};

#endif /* GEOMETRY_H_ */
