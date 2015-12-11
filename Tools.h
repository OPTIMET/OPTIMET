#ifndef TOOLS_H_
#define TOOLS_H_

#include "Spherical.h"
#include "Cartesian.h"
#include "types.h"
#include <complex>

namespace optimet {
  //! This should be done with eigen blocks...
  template<class SCALAR, class T>
  void pushToMatrix(
      SCALAR **T_small_, t_uint rows_, t_uint columns_,
      Eigen::MatrixBase<T> &T_large_,  t_uint row_index_, t_uint column_index_) {
    typedef typename T::Scalar Scalar;
    for(t_uint i=0; i < rows_; i++)
      for(t_uint j=0; j < columns_; j++)
        T_large_(i + row_index_, j + column_index_) = static_cast<Scalar>(T_small_[i][j]);
  }
}
/**
 * The Tools class implements several static tools for common use.
 */
class Tools {
public:
  /**
   * Default constructor for the Tools class.
   */
  Tools();

  /**
   * Default destructor for the Tools class.
   */
  virtual ~Tools();

  /**
   * Calculates the distance between two points in spherical coordinates.
   * @param point1 the coordinates of point1 in Spherical<double>.
   * @param point2 the coordinates of point2 in Spherical<double>.
   * @return the distance between the two points.
   */
  static double findDistance(Spherical<double> point1,
                             Spherical<double> point2);

  /**
   * Calculates the distance between two points in cartesian coordinates.
   * @param point1 the coordinates of point1 in Cartesian<double>.
   * @param point2 the coordinates of point2 in Cartesian<double>.
   * @return the distance between the two points.
   */
  static double findDistance(Cartesian<double> point1,
                             Cartesian<double> point2);

  /**
   * Projects a spherical point onto a SphericalP vector (cartesian projection).
   * @param point the point to be projected.
   * @param vector the vector to be used as basis for projection.
   * @return the projected vector in SphericalP.
   */
  static SphericalP<std::complex<double>>
  toProjection(Spherical<double> point,
               SphericalP<std::complex<double>> vector);

  /**
   * Projects a spherical point onto a SphericalP vector (spherical projection).
   * @param point the point to be projected.
   * @param vector the vector to be used as basis for projection.
   * @return the projected vector in SphericalP.
   */
  static SphericalP<std::complex<double>>
  fromProjection(Spherical<double> point,
                 SphericalP<std::complex<double>> vector);

  /**
   * Converts a point from spherical to SphericalP coordinates.
   * @param point - the coordinates of the point in
   * Spherical<std::complex<double> >.
   * @return the coordinates of the point in SphericalP<std::complex<double> >.
   */
  static SphericalP<std::complex<double>>
  toSphericalP(Spherical<std::complex<double>> point);

  /**
   * Converts a point from spherical to SphericalP coordinates.
   * @param point - the coordinates of the point in Spherical<double>.
   * @return the coordinates of the point in SphericalP<double>.
   */
  static SphericalP<double> toSphericalP(Spherical<double> point);

  /**
   * Converts a point from spherical to cartesian coordinates.
   * @param point - the coordinates of the point in Spherical<double>.
   * @return the coordinates of the point in Cartesian<double>.
   */
  static Cartesian<double> toCartesian(Spherical<double> point);

  /**
   * Converts a point from cartesian to spherical coordinates.
   * @param point - the coordinates of the point in Cartesian<double>.
   * @return the coordinates of the point in Spherical<double>.
   */
  static Spherical<double> toSpherical(Cartesian<double> point);

  /**
   * Compares two doubles to see if they are "equal".
   * The internal parameter maxUlps specifies the maximum number of doubles
   * that can be generated between A and B before they are considered
   * non-equal.
   * @param A the first double.
   * @param B the second double.
   * @return true if equal, false if not.
   */
  static bool equalDoubles(double A, double B);

  /**
   * Returns a compound iterator from two m and n iterators.
   * @param m the first iterator.
   * @param n the second iterator.
   * @return the compound iterator.
   */
  static long toCompound(long m, long n);

  /**
   * Returns the first single iterator (m) from a compound one.
   * @param p the compound iterator.
   * @return the first iterator.
   */
  static long toFirstIt(long p);

  /**
   * Returns the second single iterator (n) from a compound one.
   * @param p the compound iterator
   * @return the second iterator.
   */
  static long toSecondIt(long p);

  /**
   * Calculates the maximum value of a compound iterator with n.
   * @param n the maximum value of n.
   * @return the maximum value of a compound iterator with n.
   */
  static long iteratorMax(long n);

  /**
   * Build a unit matrix.
   * @param size_ the size of the matrix (rows or columns).
   * @param I_ the unit matrix.
   */
  static void makeUnitMatrix(long size_, std::complex<double> **I_);

  /**
   * Push the T_small_ matrix into the T_large_ matrix.
   * @param T_small_ the matrix to be pushed.
   * @param rows_ the number of rows of T_small_.
   * @param columns_ the number of columns of T_small_.
   * @param T_large_ the host matrix.
   * @param row_index_ the row index to push in.
   * @param column_index_ the column index to push in.
   */
  static void pushToMatrix(std::complex<double> **T_small_, long rows_,
                           long columns_, std::complex<double> **T_large_,
                           long row_index_, long column_index_);

  /**
   * Returns the value of an Associated Legendre function \f$P_n^m(x)\f$.
   * Static member. Insures m < 0 compatibility.
   * @param argument_ the argument x.
   * @param orderN_ the order n.
   * @param orderM_ the order m.
   * @return the Associated Legendre function.
   */
  static double getLegendre(double argument_, int orderN_, int orderM_);

  /**
   * Allocate space for a 2D matrix.
   * @param i_x the 1st size of the matrix.
   * @param i_y the 2nd size of the matrix.
   * @return pointer to the allocated matrix.
   */
  static std::complex<double> **Get_2D_c_double(int i_x, int i_y);

  /**
   * Allocate space for a 4D matrix.
   * @param t the 1st size of the matrix.
   * @param x the 2nd size of the matrix.
   * @param y the 3rd size of the matrix.
   * @param z the 4th size of the matrix.
   * @return pointer to the allocated matrix.
   */
  static std::complex<double> ****Get_4D_c_double(int t, int x, int y, int z);

  /**
   * Translate from one coordinate system to another.
   * @param R the point to be translated
   * @param P the center of the new coordinate system.
   * @return the new coordinate of point R in center P.
   */
  static Spherical<double> toPoint(Spherical<double> R, Spherical<double> P);

  /**
   * Convert polar coordinates to cartesian coordinates.
   * @param rho_ the radial polar coordinate.
   * @param theta_ the angular polar coordinate.
   * @param x_ the x cartesian coordinate.
   * @param y_ the y cartesian coordinate.
   */
  static void Pol2Cart(double rho_, double theta_, double &x_, double &y_);
};

#endif /* TOOLS_H_ */
