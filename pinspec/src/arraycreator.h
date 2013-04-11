/**
 * @file arraycreator.h
 * @brief Functions for creating arrays with similar syntax to MATLAB
 * @author William Boyd (wboyd@mit.edu)
 * @date March 28, 2012 
 *
 */

#ifndef ARRAYCREATOR_H_
#define ARRAYCREATOR_H_

#ifdef __cplusplus
#include <limits>
#endif

/**
 * @brief Creates an array of equally spaced values
 * @details Helper function to generate an array of equally spaced values 
 *          between a specified start and end point. Modeled after MATLAB's 
 *          linspace function.
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
 * @brief Creates an array of equal logarithmically spaced values
 * @details Helper function to generate an array of equal logarithmically 
 *          spaced values between a specified start and end point. Modeled 
 *          after MATLAB's logspace function.
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
