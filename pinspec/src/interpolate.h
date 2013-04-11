/**
 * @file interpolate.h
 * @brief Utility functions for linear interpolation of 1D functions
 * @author William Boyd (wboyd@mit.edu)
 * @date March 13, 2012 
 *
 */

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#ifdef __cplusplus
#include <limits>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#endif

/**
 * @brief This function finds the index of the first element in an array that 
 *        is greater than the given parameter value. This is a recursive
 *        function and it uses the binary search algorithm to find the upper 
 *        bound index.
 * @param x array of values to search
 * @param upper_bound the current upper bound at this level of recursion
 * @param lower_bound the current lower bound at this level of recursion
 * @param pt the value to search for
 * @return the index of the first element in x which is greater than pt
 */
template <typename T, typename U>
int findUpperIndex(T* x, int upper_bound, int lower_bound, U pt) {

     /* Compute the delta between the two bounding indices */
     int bound_delta = upper_bound - lower_bound;

    /* Compute the midpoint between the bounding indices */
    int new_bound = floor(bound_delta / 2.0) + lower_bound;

    /* Check that bound are appropriate - if not, return infinity */
    if (bound_delta <= 0)
        return std::numeric_limits<int>::infinity();

    /* If the upper and lower bound only differ by one, we are are finished */
    if (bound_delta == 1)
        return upper_bound;

    /* If the pt is below the new bound, the new bound becomes upper bound */
    else if (pt <= x[new_bound])
        return findUpperIndex(x, new_bound, lower_bound, pt);

    /* If the pt is above the new bound, the new bound becomes lower bound */
    else
         return findUpperIndex(x, upper_bound, new_bound, pt);
}


/**
 * @brief This function takes in the x and y values of a 1D function and 
 *        returns the linearly interpolated y value at a particular x-coordinate
          (pt). It assumes that the values given in x are ordered from least to 
 *        greatest.
 * @param x vector of x values
 * @param y vector of y values
 * @param length the number of x and y
 * @param pt x-coordinate we wish to interpolate
 */

template <typename T, typename U, typename P>
P linearInterp(T* x, T* y, int length, U pt) {

    int index = 0;
    P y_pt = 0;

    /* If the length given is less than zero, exit program */
    if (length <= 0)
      exit(1);

    /* If the length is exactly 1 then return the only y value */
    else if (length == 1)
         y_pt =  y[0];

    /* If the pt is less than all x values, return the y corresponding to
     * the least x value */
    else if (pt <= x[0])
        y_pt =  y[0];

    /* If the pt is greater than all x values, return the y corresponding to
     * the greatest x value */
    else if (pt >= x[length-1])
        y_pt =  y[length-1];

    /* Otherwise the pt is sandwiched between two of the x values */
    else {
        /* Find the index of x which bounds pt */
	index = findUpperIndex(x, length-1, 0, pt);

	/* Find the slope of the line between the two sandwich points */
	double m = (y[index] - y[index-1]) / (x[index] - x[index-1]);

	/* Return the interpolated point using the point-slope formula */
	y_pt =  m * (pt - x[index]) + y[index];
    }

    return y_pt;
}


#endif /* INTERPOLATE_H_ */
