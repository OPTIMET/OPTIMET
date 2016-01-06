#ifndef SPHERICALP_H_
#define SPHERICALP_H_

/**
 * The SphericalP class implements shpherical projection coordinates.
 * This is a template class so any standard type may be used to create
 * spherical coordinates.
 * @warning Do not use without initialization.
 */
template <class carType> class SphericalP {
public:
  carType rrr; /**< The uR coordinate. */
  carType the; /**< The uTheta coordinate. */
  carType phi; /**< The uPhi coordinate. */

  /**
   * Initialization constructor for SphericalP.
   * @param rrr_ the carType value for urrr
   * @param the_ the carType value for uthe
   * @param phi_ the carType value for uphi
   * @see init()
   */
  SphericalP(carType rrr_, carType the_, carType phi_) {
    init(rrr_, the_, phi_);
  }

  /**
   * Default SphericalP constructor.
   */
  SphericalP(void) {
    //
  }

  /**
   * Default SphericalP destructor.
   */
  ~SphericalP(void) {
    //
  }

  /**
   * Initialize a SphericalP object.
   * @param rrr_ the carType value for x
   * @param the_ the carType value for y
   * @param phi_ the carType value for z
   * @see Cartesian()
   */
  void init(carType rrr_, carType the_, carType phi_) {
    rrr = rrr_;
    the = the_;
    phi = phi_;
  }

  /**
   * Implements the dot-product of two spherical projection vectors.
   * @param argument_ the vector to be multiplied with.
   * @return the dot-product of this and argument_.
   */
  carType operator*(SphericalP<carType> argument_) {
    return rrr * argument_.rrr + the * argument_.the + phi * argument_.phi;
  }

  SphericalP<carType> operator+(SphericalP<carType> argument_) {
    return SphericalP<carType>(rrr + argument_.rrr, the + argument_.the,
                               phi + argument_.phi);
  }

  SphericalP<carType> operator-(SphericalP<carType> argument_) {
    return SphericalP<carType>(rrr - argument_.rrr, the - argument_.the,
                               phi - argument_.phi);
  }

  SphericalP<carType> operator*(carType argument_) {
    return SphericalP<carType>(rrr * argument_, the * argument_,
                               phi * argument_);
  }
};

#endif /* SPHERICALP_H_ */
