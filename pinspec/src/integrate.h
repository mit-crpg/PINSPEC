/**
 * @file integrate.h
 * @brief Utility functions for numerical integration of 1D functions.
 * @details This module provides a set of templated integration methods for both
 *          normal numerical integration and cumulative numerical integration.
 *          It implements the Riemann integration method with left-aligned, 
 *          right-aligned and midpoint methods, the trapezoidal integration, 
 *          and higher-order Simpsons, Simpson's 3/8ths and Boole's methods.
 * @author William Boyd (wboyd@mit.edu)
 * @date February 9, 2012 
 *
 */

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#ifdef __cplusplus
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#endif


/**
 * @enum integrationSchemes
 * @brief Numerical integration methods
 */

/**
 * @var integrationScheme
 * @brief Numerical integration methods
 */
typedef enum integrationSchemes {
    /** Left-aligned Riemann integration method */
    RIEMANN_LEFT,
    /** Right-aligned Riemann integration method */	
    RIEMANN_RIGHT,
    /** Midpoint Riemann integration method */
    RIEMANN_CENTER,
    /** Trapezoidal integration method */
    TRAPEZOIDAL,
    /** Simpson's integration method */
    SIMPSONS,
    /** Simpson's 3/8ths integration method */
    SIMPSONS38,
    /** Boole's integration method */
    BOOLES,
} integrationScheme;


/**
 * @brief Performs a cumulative numerical integral over arrays of x and y 
          values using the specificed integration method.
 * @param x the the x values
 * @param y the y values
 * @param cdf the array of cdf values at each value of x and y
 * @param length the length of the x and y arrays
 * @param scheme the integration scheme
 */
template <typename T, typename U>
void cumulativeIntegral(T* x, T* y, U* cdf, int length, 
			integrationScheme scheme) {

    /* Calculate cumulative integral */
    for (int i=1; i < length+1; i++)
        cdf[i-1] = (U)integrate(x, y, i, scheme);

    return;
}


/**
 * @brief Performs a numerical integral using the left-centered Riemann
 *        numerical integration method.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template <typename T>
double computeRiemannLeft(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 0; i < length; i++) {
        if (i < length - 1)
	    delta_x = x[i+1] - x[i];
	else
	    delta_x = x[i] - x[i-1];

        integral += delta_x * y[i];
    }

    return integral;
}


/**
 * @brief Performs a numerical integral using the right-centered Riemann
 *        numerical integration method.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template <typename T>
double computeRiemannRight(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 1; i < length-1; i++) {
        delta_x = x[i] - x[i-1];
	integral += delta_x * y[i];
    }

    return integral;
}


/**
 * @brief Performs a numerical integral using the midpoint-centered Riemann 
 *        numerical integration method.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template <typename T>
double computeRiemannCenter(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 0; i < length; i++) {
        if (i == 0)
	    delta_x = x[i+1] - x[i];
        else if (i == length - 1)
	    delta_x = x[i] - x[i-1];
	else
	    delta_x = (x[i] - x[i-1]) / 2.0 + (x[i+1] - x[i]) / 2.0;

	integral += delta_x * y[i];
    }

    return integral;
}


/**
 * @brief Performs a numerical integral using the trapezoidal numerical
 *        integration method.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template <typename T>
double computeTrapezoidal(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 1; i < length; i++) {
        delta_x = x[i] - x[i-1];
	integral += delta_x * (y[i] + y[i-1]) / 2.0;
    }

    return integral;
}


/**
 * @brief Performs a numerical integral using the Simpson's numerical 
 *        integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template <typename T>
double computeSimpsons(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 0; i < length - 2; i++) {
        delta_x = x[i+1] - x[i];
	integral += (delta_x / 6.0) * (y[i] + 4*y[i+1] + y[i+2]);
    }

    return integral;
}


/**
 * @brief Performs a numerical integral using the Simpson's 3/8ths
 *        numerical integration method.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template <typename T>
double computeSimpsons38(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 0; i < length - 3; i++) {
        delta_x = x[i+1] - x[i];
	integral += (delta_x / 8.0) * (y[i] + 3*y[i+1] + 3*y[i+2] + y[i+3]);
    }

    return integral;
}


/**
 * @brief Performs a numerical integral using the Boole's numerical
 *        integration method.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @return the integral value
 */
template<typename T>
double computeBooles(T* x, T* y, int length) {

    double integral = 0;
    double delta_x = 0;

    for (int i = 0; i < length - 4; i++) {
        delta_x = x[i+1] - x[i];
	integral += (delta_x / 90.0) * (7*y[i] + 32*y[i+1] + 12*y[i+2] +
					          32*y[i+3] + 7*y[i+4]);
    }

    return integral;
}


/**
 * @brief Performs a 1D numerical integral over arrays of x and y values using 
          the specified integration scheme.
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @param scheme integration scheme
 * @return the integral value
 */
template <typename T>
double integrate(T* x, T* y, int length, integrationScheme scheme) {

    double integral = 0;

    switch(scheme) {
        case RIEMANN_RIGHT:
	    integral = computeRiemannRight(x, y, length);
	    break;
	case RIEMANN_LEFT:
	    integral = computeRiemannLeft(x, y, length);
	    break;
        case RIEMANN_CENTER:
	    integral = computeRiemannCenter(x, y, length);
	    break;
	case TRAPEZOIDAL:
	    integral = computeTrapezoidal(x, y, length);
	    break;
	case SIMPSONS:
	    integral = computeSimpsons(x, y, length);
	    break;
	case SIMPSONS38:
	    integral = computeSimpsons38(x, y, length);
	    break;
	case BOOLES:
	    integral = computeBooles(x, y, length);
	    break;
	}

    return integral;
}


#endif
