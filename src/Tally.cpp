/*
 * tally.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Tally.h"


/**
 * Default Tally constructor
 */
Tally::Tally() {

	_name = (char*)"";

	 /* Sets the default delta between bins to zero */
	_bin_delta = 0;

	/* Default is to tally all isotopes */
	_isotopes = (char*)"all";
}


/**
 * Tally destructor deletes memory for tallies, number of tallies,
 * bin centers and bin edges if they have been created
 */
Tally::~Tally() {

	if (_num_bins != 0) {
		delete [] _tallies;
		delete [] _num_tallies;
		delete [] _centers;
		if (_bin_spacing != OTHER)
			delete [] _edges;
	}
}


/**
 * Returns the name of this Tally as specified by the user
 * @return the Tally's name
 */
char* Tally::getTallyName() {
	return _name;
}


/**
 * Returns the number of bins
 * @return the number of bins
 */
int Tally::getNumBins() {
	return _num_bins;
}


/**
 * Returns a double array of bin edge values
 * @return array of bin edge values
 */
float* Tally::getBinEdges() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin edges for Tally %s since "
				 "the bins have not yet been created", _name);

	 return _edges;
}


/**
 * Returns a double array of bin center values
 * @return array of bin center values
 */
double* Tally::getBinCenters() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin centers for Tally %s since the "
				 "bins have not yet been created", _name);

	 return _centers;
}


/**
 * Returns the delta spacing between bins. NOTE: this value is only non-zero
 * for EQUAL and LOGARITHMIC bin types
 * @return the spacing between bins
 */
float Tally::getBinDelta() {
	return _bin_delta;
}


float Tally::getBinDelta(float sample) {

	/* If this Tally uses equally spaced bins in linear or logarithmic
	 * space, return the bin delta */
	if (_bin_spacing == EQUAL || _bin_spacing == LOGARITHMIC)
		return _bin_delta;

	/* If instead this Tally uses irregularly spaced bin edges defined
	 * by a user, compute bin delta of the bin around the sample */
	else {
		int bin_index = getBinIndex(sample);
		return (_edges[bin_index] - _edges[bin_index-1]);
	}
}


/**
 * Returns the bin spacing type (EQUAL, LOGARITHMIC, OTHER)
 * @return the bin spacing type
 */
binSpacingType Tally::getBinSpacingType() {
	return _bin_spacing;
}


/**
 * Returns the type of tally for these bins (ISOTOPE, MATERIAL, REGION)
 * @return the tally type
 */
tallyDomainType Tally::getTallyDomainType() {
	return _tally_domain;
}

/**
 * Returns the type of tally for these bins (FLUX, COLLISION, ABSORPTION, etc)
 * @return the tally type
 */
tallyType Tally::getTallyType() {
	return _tally_type;
}


/**
 * Returns a double array of the tallies within each bin
 * @return an array of
 */
double* Tally::getTallies() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tallies for Tally %s since the "
				 "bins have not yet been created", _name);

	 return _tallies;
}


/**
 * Returns a specific tally for a specific bin
 * @param bin_index the index for the bin of interest
 * @return the tally within that bin
 */
double Tally::getTally(int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally for a bin index for Tally %s"
				"which does not exist: %d, num_bins = %d", _name, bin_index,
				_num_bins);

	return _tallies[bin_index];
}


/**
 * Returns an int array of the number of times tallied within each bin
 * @return an array of the number of tallies in each bin
 */
int* Tally::getNumTallies() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tally numbers for Tally %s since "
				 "the bins have not yet been created", _name);

	return _num_tallies;
}


/**
 * Returns the number of times tallied within a specific bin
 * @param bin_index the bin of interest
 * @return the number of tallies in that bin
 */
int Tally::getNumTallies(int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally number for Tally %s for "
				"a bin index which does not exist: %d", _name, bin_index);

	return _num_tallies[bin_index];
}


/**
 * Returns the maximum tally value over all bins
 * @return the maximum tally value
 */
double Tally::getMaxTally() {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return the maximum tally for Tally %s"
				 "since the bins have not yet been created", _name);

	double max_tally = 0;

	/* Loop over all bins */
	for (int i=0; i < _num_bins; i++) {
		if (_tallies[i] > max_tally)
			max_tally = _tallies[i];
	}

	return max_tally;
}


/**
 * Returns the maximum tally value over all bins
 * @return the maximum tally value
 */
double Tally::getMinTally() {
	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return the minimum tally for Tally %s"
				 " since the bins have not yet been created", _name);

	double min_tally = std::numeric_limits<double>::infinity();

	/* Loop over all bins */
	for (int i=0; i < _num_bins; i++) {
		if (_tallies[i] < min_tally)
			min_tally = _tallies[i];
	}

	return min_tally;
}


/**
 * Finds the bin index for a sample in a set of bins. If the samples
 * is outside the bounds of all bins, it returns infinity
 * @param sample the sample value of interest
 * @return the bin index for the sample
 */
int Tally::getBinIndex(float sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return a bin index for Tally %s since "
				 "the bins have not yet been created", _name);

	/* Set index to infinity to begin with */
	int index = std::numeric_limits<float>::infinity();

	/* if the sample is equal to the last bin edge, return the last bin */
	if (sample == _edges[_num_bins])
		return _num_bins-1;

	/* Equally spaced bins */
	if (_bin_spacing == EQUAL)
		index = int((sample - _edges[0]) / _bin_delta);

	/* Logarithmically spaced bins */
	else if (_bin_spacing == LOGARITHMIC)
		index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

	/* If the bin_type == OTHER then the bin edges were not generated by
	 * generateEqualBinEdges, so use a brute force search to find the bin */
	else {

		/* Loop over all bin edges to find the correct bin index */
		for (int i=0; i <= _num_bins; i++) {
			if (sample >= _edges[i] && sample < _edges[i+1]) {
				index = i;
				break;
			}
		}
	}

	/* If this sample was not contained within a bin set index to infinity*/
	if (index > _num_bins)
		index = std::numeric_limits<float>::infinity();

	return index;
}


/**
 * Return the isotopes that this Tally is meant to tally
 * @return a character array of the isotope's name or "all" for all isotopes
 */
char* Tally::getIsotopes() {
	return _isotopes;
}



/**
 * Sets this Tally's name
 * @param name the name of the Tally
 */
void Tally::setTallyName(char* name) {
	_name = name;
}


/**
 * Set the domain type for this Tally (ISOTOPE, MATERIAL, REGION)
 * @param type the tally type
 */
void Tally::setTallyDomainType(tallyDomainType type) {
	_tally_domain = type;
}


/**
 * Set the type of tally for this Tally (FLUX, COLLISION, ABSORPTION)
 * @param type the tally type
 */
void Tally::setTallyType(tallyType type) {
	_tally_type = type;
}


/**
 * Set a user-defined double array of bin edge values
 * @param edges the array of bin edges
 * @param num_bins the number of bins
 */
void Tally::setBinEdges(float* edges, int num_bins) {

	_num_bins = num_bins;
	_edges = edges;
	_bin_spacing = OTHER;

	/* Set all tallies to zero by default */
	_tallies = new double[num_bins];
	_num_tallies = new int[num_bins];


	/* Loop over tallies and set to zero */
	for (int i=0; i < _num_bins; i++) {
		_tallies[i] = 0.0;
		_num_tallies[i] = 0;
	}

	/* Create an array of the center values between bins */
	generateBinCenters();
}



/**
 * Set the isotope that this Tally is meant to tally
 * @param isotopes a character array of the isotope's name or
 * "all" for all isotopes
 */
void Tally::setIsotopes(char* isotopes) {
	_isotopes = isotopes;
}


/**
 * Generate edges between bins defined by a start and end point
 * @param start first bin edge value
 * @param end last bin edge value
 * @param num_bins the number of bins to be created
 * @param type the type of bins (EQUAL or LOGARITHMIC)
 */
void Tally::generateBinEdges(float start, float end, int num_bins,
												binSpacingType type) {
	if (start == end)
		log_printf(ERROR, "Unable to create bins for Tally %s between"
				"the same start and end points: %f", _name, start);

	_num_bins = num_bins;
	_bin_spacing = type;

	/* Allocate memory for tallies */
	_tallies = new double[num_bins];
	_num_tallies = new int[num_bins];

	/* Set all tallies to zero by default */
	for (int i=0; i < num_bins; i++) {
		_tallies[i] = 0;
		_num_tallies[i] = 0;
	}

	/* Equal spacing between bins */
	if (type == EQUAL) {
		_bin_delta = float(end - start) / float(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = linspace<float, float>(start, end, num_bins+1);
	}

	/* Logarithmically equal spacing between bins */
	else if (type == LOGARITHMIC) {
		_bin_delta = float(log10(end) - log10(start)) / float(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = logspace<float, float>(start, end, num_bins+1);
	}

	else
		log_printf(ERROR, "Bin type %d is not yet implemented for Tally %s",
															_name, type);

	/* Create an array of the center values between bins */
	generateBinCenters();

	return;
}


/**
 * Compute the center points between bin edges for this Tally's bins
 */
void Tally::generateBinCenters() {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot generate bin centers for Tally %s since "
				 "the bins have not yet been created", _name);

	/* Allocate memory for the bin centers array */
	_centers = new double[_num_bins];

	/* Loop over all bins and find the midpoint between edges */
	for (int i=0; i < _num_bins; i++)
		_centers[i] = (_edges[i] + _edges[i+1]) / 2.0;

	return;
}


/**
 * Tallies unity for each sample in a double array of samples
 * @param samples array of samples to tally
 * @param num_samples the number of samples to tally
 */
void Tally::tally(float* samples, int num_samples) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally samples in Tally %s since the "
				 "bins have not yet been created", _name);

	int bin_index;

	/* Loop over and tally all samples */
	for (int i=0; i < num_samples; i++) {
		bin_index = getBinIndex(samples[i]);
		if (bin_index >= 0 && bin_index < _num_bins) {
			_tallies[bin_index]++;
			_num_tallies[bin_index]++;
		}
	}

	return;
}


/**
 * Tallies unity for a sample
 * @param samples array of samples to tally
 */
void Tally::tally(float sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally sample in Tally %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(sample);

	if (bin_index >= 0 && bin_index < _num_bins) {
		_tallies[bin_index]++;
		_num_tallies[bin_index]++;
	}

	return;
}


/**
 * Tallies a weight for each sample in a double array of samples
 * @param samples array of samples to tally
 * @param sample_weights array of sample weights to increment tallies by
 * @param num_samples the number of samples to tally
 */
void Tally::weightedTally(float* samples, float* sample_weights,
														int num_samples) {
	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted samples in Tally %s "
				 "since the bins have not yet been created", _name);

	int bin_index;

	/* Loop over and tally all samples */
	for (int i=0; i < num_samples; i++) {
		bin_index = getBinIndex(samples[i]);
		if (bin_index >= 0 && bin_index < _num_bins) {
			_tallies[bin_index] += sample_weights[i];
			_num_tallies[bin_index]++;
		}
	}

	return;
}


/**
 * Tallies a weight for a sample
 * @param sample a sample to tally
 * @param weight the weight to increment tally by
 */
void Tally::weightedTally(float sample, float weight) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Tally %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(sample);

	if (bin_index >= 0 && bin_index < _num_bins) {
		_tallies[bin_index] += double(weight);
		_num_tallies[bin_index]++;
	}

	return;
}


/**
 * Divide each tally by the maximum tally value
 */
void Tally::normalizeTallies() {

	if (_num_bins == 0)
		log_printf(ERROR, "Cannot normalize tallies for Tally %s since it is"
						"the bins have not yet been created", _name);

	double max_tally = getMaxTally();

	/* Divide each tally by maximum tally value */
	for (int n=0; n < _num_bins; n++)
		_tallies[n] /= max_tally;

	return;
}


/**
 * Divide each tally by a given scaling factor
 * @param scale_factor factor to normalize tallies by
 */
void Tally::normalizeTallies(float scale_factor) {

	if (_num_bins == 0)
		log_printf(ERROR, "Cannot normalize tallies for Tally %s since it is"
						"the bins have not yet been created", _name);

	/* Divide each tally by maximum tally value */
	for (int n=0; n < _num_bins; n++)
		_tallies[n] /= double(scale_factor);

	return;
}
