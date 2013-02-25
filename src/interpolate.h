/*
 * interpolate.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


/**
 * This function takes in the x and y values of a function and returns the
 * linearly interpolated y values at an particular x (pt). It assumes
 * that the values given in x are ordered from least to greatest
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
	if (length <= 0) {
		printf("Tried to interpolate on a vector of length 0\n");
		exit(1);
	}

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


/**
 * This function takes in the x and y values of a function and returns the
 * cubic spline interpolated y values at an particular x (pt). It assumes
 * that the values given in x are ordered from least to greatest
 * @param x vector of x values
 * @param y vector of y values
 * @param length the number of x and y
 * @param pt x-coordinate we wish to interpolate
 */
template <typename T, typename U, typename P>
P splineInterp(T* x, T* y, int length, U pt) {
	P y_pt = 0;

	//TODO
	printf("Spline interpolation is not yet implemented\n");
	exit(1);

	return y_pt;
}


/**
 * This function finds the index of the first element in an array that is
 * greater than the given parameter value. This is a recursive function and
 * it uses the binary search algorithm to find the upper bound index.
 * @param x double array we wish to search
 * @param upper_bound the current upper bound at this level of recursion
 * param lower_bound the current loewr bound at this level of recursion
 * @param pt the point of interest
 * @return the index of the first element in x which is greater than pt
 */
template <typename T, typename U>
int findUpperIndex(T* x, int upper_bound, int lower_bound, U pt) {

	/* Compute the delta between the two bounding indices */
	int bound_delta = upper_bound - lower_bound;

	/* Compute the midpoint between the bounding indices */
	int new_bound = floor(bound_delta / 2.0) + lower_bound;

	/* Check that bound are appropriate */
	if (bound_delta <= 0) {
		printf("Unable to find the upper index using binary search"
				"since upper_bound = %d and lower_bound = %d",
				upper_bound, lower_bound);
		exit(1);
	}

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


#endif /* INTERPOLATE_H_ */
