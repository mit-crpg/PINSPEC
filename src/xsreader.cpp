/*
 * xsreader.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "xsreader.h"


/**
 * This method parses and input file of cross-section data and loads the
 * energy values (as eV) and the cross-section values (barns) into two
 * float arrays. It takes in the energies in MeV since that is the
 * most typical unit for ENDF nuclear data, but converts the values
 * into eV during parsing.
 * @param file the filename for the data
 * @param energies a pointer to a float array for the energies (MeV)
 * @param xs_values a pointer to a float array fo the xs values (barns)
 * @param delimiter the type of delimiter between energy and xs
 * @return the number of data points
 */
int parseCrossSections(const char* file, float* energies, float* xs_values) {

	/* Instantiate I/O variables */
	std::ifstream input_file(file, std::ios::in);
	std::string buff;
	int count = 0;

	/* get the first line of the input file */
	getline(input_file, buff);

	/* Parse over each line in the file */
	while(getline(input_file, buff)) {
		/* Load data into the arrays for this line */
		sscanf(buff.c_str(), "%f %f", &energies[count], &xs_values[count]);
		count++;
	}

	/* Convert energy values from MeV to eV */
	for (int i=0; i < count; i++)
		energies[i] *= 1E6;

	/* Close the file and return the number of data points */
	input_file.close();

	return count;
}


/**
 * This function counts the number of lines in an input file with
 * nuclear cross-sections
 * @param filename the file of interest
 * @return the number of lines in the file
 */
int getNumCrossSectionDataPoints(const char* filename) {

	/* Instantiate I/O variables */
	std::ifstream input_file(filename, std::ios::in);
	std::string line;
	int num_xs_values = 0;

	/* Loop over each line and update counter */
	while(getline(input_file, line))
	    num_xs_values++;

	/* Close the file and return the number of lines */
	input_file.close();
	return num_xs_values;
}

