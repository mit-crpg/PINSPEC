/**
 * @file vector.h
 * @brief Utility functions for commonly used vector operations.
 * @author William Boyd (wboyd@mit.edu)
 * @date April 14, 2013
 *
 */

#ifndef VECTOR_H_
#define VECTOR_H_


/**
 * @brief Computes the dot product of two vectors in the 2D cartesian space.
 * @details This method computes the scalar product of two 2D vectors as follows:
 *       \f$ \vec{u}\cdot\vec{v}=\vec{u}_x *\vec{v}_x + \vec{u}_y * \vec{u}_y \f$
 * @param x1 the x-coodinate of the first vector
 * @param y1 the y-coordinate of the first vector
 * @param x2 the x-coordinate of the second vector
 * @param y2 the y-coordinate of the second vector
 */
template <typename T>
T dotProduct2D(T x1, T y1, T x2, T y2) {
    return x1 * x2 + y1 * y2;
}


/**
 * @brief Computes the dot product of two vectors in 3D cartesian space.
 * @details This method computes the scalar product of two 3D vectors as follows:
 *       \f$ \vec{u}\cdot\vec{v}=\vec{u}_x\cdot\vec{v}_x + 
 *       \vec{u}_y\cdot\vec{u}_y + \vec{u}_z\cdot\vec{v}_z \f$
 * @param x1 the x-coodinate of the first vector
 * @param y1 the y-coordinate of the first vector
 * @param z1 the z-coordinate of the first vector
 * @param x2 the x-coordinate of the second vector
 * @param y2 the y-coordinate of the second vector
 * @param z2 the z-coordinate of the second vector
 */
template <typename T>
T dotProduct3D(T x1, T y1, T z1, T x2, T y2, T z2) {
    return x1 * x2 + y1 * y2 + z1 * z2;
}


/**
 * @brief Computes the norm of a vector in the 2D cartesian space.
 * @details This method computes the norm of a vector \f$ \vec{u} \f$ in 2D 
 *          as follows:
 *          \f$ \| \vec{u} \| = \sqrt{x \cdot x + y \cdot y} \f$
 * @param x the x-coordinate of the vector
 * @param y the y-coordinate of the vector
 */
template <typename T>
T norm2D(T x, T y) {
  return sqrt(x * x + y * y);
}


/**
 * @brief Computes the norm of a vector in 3D cartesiance space.
 * @details This method computes the norm of a vector \f$ \vec{u} \f$ in 3D 
 *          as follows:
 *          \f$ \| \vec{u} \| = \sqrt{x\ cdot x + y \cdot y + z \cdot z} \f$
 * @param x the x-coordinate of the vector
 * @param y the y-coordinate of the vector
 * @param z the z-coordinate of the vector
 */
template <typename T>
T norm3D(T x, T y, T z) {
  return sqrt(x * x + y * y + z * z);
}


#endif /* VECTOR_H_ */
