/*
 * arraycreator.h
 *
 *  Created on: Mar 28, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef ARRAYCREATOR_H_
#define ARRAYCREATOR_H_

/* A set of templated helper functions inspired by MATLAB */

/**
 * Helper function to generate an array of equally spaced doubles between
 * a given start and end point. Modeled after MATLAB's linspace function
 * @param start the starting point
 * @param end the ending point
 * @param num_values the number of values to create
 * @return a pointer to the array of points
 */
template <typename T, typename U>
U* linspace(T start, T end, int num_values) {

	U* values = new U[num_values];

	/* Spacing between values */
	double delta = double(end - start) / double(num_values-1);

	/* Loop over all values */
	for (int i=0; i < num_values; i++)
		values[i] = delta * i + start;

	return values;
}


/**
 * Helper function to generate an array of equal logarithmically spaced
 * doubles between a given start and end point. Modeled after MATLAB's
 * linspace function
 * @param start the starting point
 * @param end the ending point
 * @param num_values the number of values to create
 * @return a pointer to the array of points
 */
template <typename T, typename U>
U* logspace(T start, T end, int num_values) {

	/* Create an equally spaced array of base 10 exponent values */
	U* values = linspace<T,U>(log10(start), log10(end), num_values);

	/* Loop over all values and project back original space */
	for (int i=0; i < num_values; i++)
		values[i] = pow(10, values[i]);

	return values;
}

#endif /* ARRAYCREATOR_H_ */
