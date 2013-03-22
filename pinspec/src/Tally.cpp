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
Tally::Tally(char* tally_name) {

	_tally_name = tally_name;
    _trigger_type = NONE;

	 /* Sets the default delta between bins to zero */
	_bin_delta = 0;

	/* Sets the default for batch statistics */
	_num_batches = 0;
    _num_bins = 0;
	_computed_statistics = false;
}


/**
 * Tally destructor deletes memory for tallies, number of tallies,
 * bin centers and bin edges if they have been created
 */
Tally::~Tally() {

	if (_num_bins != 0)
        delete [] _edges;

	if (_num_batches != 0) {
		delete [] _tallies;
		delete [] _centers;
		delete [] _batch_mu;
		delete [] _batch_variance;
		delete [] _batch_std_dev;
		delete [] _batch_rel_err;
	}
}


/**
 * Returns the of this Tally as specified by the user
 * @return the Tally's name
 */
char* Tally::getTallyName() {
	return _tally_name;
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
double* Tally::getBinEdges() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin edges for Tally %s since "
				 "the bins have not yet been created", _tally_name);

	 return _edges;
}


/**
 * Returns a double array of bin center values
 * @return array of bin center values
 */
double* Tally::getBinCenters() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin centers for Tally %s since the "
				 "centers have not yet been created", _tally_name);

	 return _centers;
}


/**
 * Returns the delta spacing between bins. NOTE: this value is only non-zero
 * for EQUAL and LOGARITHMIC bin types
 * @return the spacing between bins
 */
double Tally::getBinDelta() {
	return _bin_delta;
}


double Tally::getBinDelta(double sample) {

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
 * Returns the type of tally for these bins (FLUX, COLLISION_RATE, etc)
 * @return the tally type
 */
tallyType Tally::getTallyType() {
	return _tally_type;
}


/**
 * Returns a double array of the tallies within each bin
 * @return an array of
 */
double** Tally::getTallies() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tallies for Tally %s since the "
				 "bins have not yet been created", _tally_name);

	 return _tallies;
}


/**
 * Returns a specific tally for a specific bin
 * @param bin_index the index for the bin of interest
 * @return the tally within that bin
 */
double Tally::getTally(int batch_num, int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally for a bin index for Tally %s"
				"which does not exist: %d, num_bins = %d", 
						_tally_name, bin_index, _num_bins);
	if (batch_num < 0 || batch_num >= _num_batches)
		log_printf(ERROR, "Tried to get a tally for a batch for Tally %s"
				"which does not exist: %d, num_batches = %d", 
						_tally_name, batch_num, _num_batches);

	return _tallies[batch_num][bin_index];
}


/**
 * Returns the maximum tally value over all bins
 * @return the maximum tally value
 */
double Tally::getMaxTally() {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return the maximum tally for Tally %s"
				 "since the bins have not yet been created", _tally_name);

	double max_tally = 0;

	/* Loop over all bins */
	for (int i=0; i < _num_batches; i++) {
		for (int j=0; j < _num_bins; j++) {
			if (_tallies[i][j] > max_tally)
				max_tally = _tallies[i][j];
		}
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
				 " since the bins have not yet been created", _tally_name);

	double min_tally = std::numeric_limits<double>::infinity();

	/* Loop over all bins */
	for (int i=0; i < _num_batches; i++) {
		for (int j=0; j < _num_bins; j++) {
			if (_tallies[i][j] < min_tally)
				min_tally = _tallies[i][j];
		}
	}

	return min_tally;
}


/**
 * Returns the number of batches for this tally
 * @return the number of batches
 */
int Tally::getNumBatches() {
	return _num_batches;
}


/**
 * Returns a pointer to an array of batch averages if they have been
 * computed
 * @return a double array of batch averages for each bin
 */
double* Tally::getBatchMu() {

	if (!_computed_statistics)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"Tally %s so batch mu cannot be returned", _tally_name);

	return _batch_mu;
}


/**
 * Returns a pointer to an array of batch variances if they have been
 * computed
 * @return a double array of batch variances for each bin
 */
double* Tally::getBatchVariance() {

	if (!_computed_statistics)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"Tally %s so batch variance cannot be returned", _tally_name);

	return _batch_variance;
}


/**
 * Returns a pointer to an array of batch standard deviations if they have
 * been computed
 * @return a double array of batch standard deviations for each bin
 */
double* Tally::getBatchStdDev() {

	if (!_computed_statistics)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"Tally %s so batch std dev cannot be returned", _tally_name);

	return _batch_std_dev;
}


/**
 * Returns a pointer to an array of batch relative errors if they have been
 * computed
 * @return a double array of batch relative errors for each bin
 */
double* Tally::getBatchRelativeError() {

	if (!_computed_statistics)
		log_printf(ERROR, "Statistics have not yet been computed for "
		"Tally %s so batch relative error cannot be returned", _tally_name);

	return _batch_rel_err;
}


double Tally::getMaxVariance() {

    double max_variance = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_variance[i] > max_variance)
            max_variance = _batch_variance[i];
    }

    return max_variance;
}


double Tally::getMaxRelErr() {

    double max_rel_err = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_rel_err[i] > max_rel_err)
            max_rel_err = _batch_rel_err[i];
    }

    return max_rel_err;
}


double Tally::getMaxStdDev() {

    double max_std_dev = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_std_dev[i] > max_std_dev)
            max_std_dev = _batch_std_dev[i];
    }

    return max_std_dev;
}


bool Tally::isPrecisionTriggered() {

    if (_trigger_type == NONE)
        return false;
    else {

        if (_trigger_type == VARIANCE) {
            if (getMaxVariance() < _trigger_precision) {
                _trigger_type = NONE;
                return false;
            }
            else {
                log_printf(INFO, "Tally %s triggered ("
                          "variance < %1.1E) with a current variance of %1.1E",
                          _tally_name, _trigger_precision, getMaxVariance());
                return true;
            }
        }
        else if (_trigger_type == STANDARD_DEVIATION) {
            if (getMaxStdDev() < _trigger_precision) {
                _trigger_type = NONE;
                return false;
            }
            else {
                log_printf(INFO, "Tally %s triggered ("
                            "std. dev. < %1.1E) with max std. dev. = %1.1E",
                            _tally_name, _trigger_precision, getMaxStdDev());       
                return true;
            }
        }
        else if (_trigger_type == RELATIVE_ERROR) {
            if (getMaxRelErr() < _trigger_precision) {
                _trigger_type = NONE;
                return false;
            }
            else {
                log_printf(INFO, "Tally %s triggered ("
                            "rel. err. < %1.1E) with max rel. err. = %1.1E",
                            _tally_name, _trigger_precision, getMaxRelErr());   
                return true;
            }
        }
    }

    return true;
}



void Tally::retrieveTallyEdges(double* data, int num_bins) {

    /* Load all tally bin centers into array */
    for (int i=0; i < _num_bins+1; i++)
        data[i] = _edges[i];

}


void Tally::retrieveTallyCenters(double* data, int num_bins){

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve bin centers for Tally %s since"
                      " it has not yet computed batch statistics", _tally_name);

    /* Load all tally bin centers into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _centers[i];

}


void Tally::retrieveTallyMu(double* data, int num_bins) {

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve tally mu for Tally %s since"
                    " it has not yet computed batch statistics", _tally_name);
    if (_num_batches == 0)
        log_printf(ERROR, "Unable to retrieve tally mu for Tally %s since it"
              " does not know how many batches it should tally", _tally_name);

    /* Load all tally mu into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _batch_mu[i];
}

void Tally::retrieveTallyVariance(double* data, int num_bins) {

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve tally variances for Tally %s "
						"since it has not yet computed batch statistics", 
																_tally_name);
    if (_num_batches == 0)
        log_printf(ERROR, "Unable to retrieve tally variances for Tally %s "
						"since it does not know how many batches it should"
						" tally", _tally_name);

    /* Load all tally variances into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _batch_variance[i];
}


void Tally::retrieveTallyStdDev(double* data, int num_bins) {

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve tally std. dev. for Tally %s "
			"since it has not yet computed batch statistics", _tally_name);
    if (_num_batches == 0)
        log_printf(ERROR, "Unable to retrieve tally std. dev. for Tally %s "
			"since it does not know how many batches it should tally", 
																_tally_name);

    /* Load all tally standard deviations into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _batch_std_dev[i];
}


void Tally::retrieveTallyRelErr(double* data, int num_bins) {

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve tally rel. err. for Tally %s "
			"since it has not yet computed batch statistics", _tally_name);
    if (_num_batches == 0)
        log_printf(ERROR, "Unable to retrieve tally rel. err. for Tally %s "	
			"since it does not know how many batches it should tally", 
																_tally_name);

    /* Load all tally variances into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _batch_rel_err[i];
}


/**
 * Set the bin spacing type for this Tally (EQUAL, LOGARITHMIC, OTHER)
 * @param type the bin spacing type
 */
void Tally::setBinSpacingType(binSpacingType type) {
    _bin_spacing = type;
}


/**
 * Set a user-defined double array of bin edge values
 * @param edges the array of bin edges
 * @param num_bins the number of bins
 */
void Tally::setBinEdges(double* edges, int num_edges) {

	_num_bins = num_edges-1;
	_bin_spacing = OTHER;
    _edges = (double*)malloc(sizeof(double)*num_edges);

    for (int i=0; i < num_edges; i++)
        _edges[i] = edges[i];

	/* Create an array of the center values between bins */
	generateBinCenters();

    return;
}


void Tally::setPrecisionTrigger(triggerType trigger_type, float precision) {

	if (precision < 0.0)
		log_printf(ERROR, "Unable to set a negative trigger precision of %f"
					" for tally %s", precision, _tally_name);

    _trigger_type = trigger_type;
    _trigger_precision = precision;
    return;
}


/**
 * Set the number of batches for this Tally. This method also
 * allocates memory for the tallies and batch statistics arrays
 * @param num_batches the number of batches
 */
void Tally::setNumBatches(int num_batches) {
	_num_batches = num_batches;

	/* Set all tallies to zero by default */
	_tallies = (double**) malloc(sizeof(double*) * _num_batches);
	for (int i=0; i < _num_batches; i++)
		_tallies[i] = new double[_num_bins];

	/* Loop over tallies and set to zero */
	for (int i=0; i < _num_batches; i++) {
		for (int j=0; j < _num_bins; j++)
			_tallies[i][j] = 0.0;
	}

	/* Allocate memory for batch-based statistical counters */
	_batch_mu = new double[_num_bins];
	_batch_variance = new double[_num_bins];
	_batch_std_dev = new double[_num_bins];
	_batch_rel_err = new double[_num_bins];
}


void Tally::incrementNumBatches(int num_batches) {

	_num_batches += num_batches;
    double** new_tallies = (double**)malloc(sizeof(double*) * _num_batches);

    for (int i=0; i < _num_batches; i++) {

        if (i < _num_batches-num_batches)
            new_tallies[i] = _tallies[i];

        else {
            new_tallies[i] = new double[_num_bins];
	        /* Loop over tallies and set to zero */
	        for (int j=0; j < _num_bins; j++)
		        new_tallies[i][j] = 0.0;
        }
    }
    
    _tallies = new_tallies;
}


/**
 * Generate edges between bins defined by a start and end point
 * @param start first bin edge value
 * @param end last bin edge value
 * @param num_bins the number of bins to be created
 * @param type the type of bins (EQUAL or LOGARITHMIC)
 */
void Tally::generateBinEdges(double start, double end, int num_bins,
												binSpacingType type) {
	if (start == end)
		log_printf(ERROR, "Unable to create bins for Tally %s between"
				"the same start and end points: %f", _tally_name, start);

	_num_bins = num_bins;
	_bin_spacing = type;

	/* Equal spacing between bins */
	if (type == EQUAL) {
		_bin_delta = double(end - start) / double(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = linspace<double, double>(start, end, num_bins+1);
	}

	/* Logarithmically equal spacing between bins */
	else if (type == LOGARITHMIC) {
		_bin_delta = double(log10(end) - log10(start)) / double(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = logspace<double, double>(start, end, num_bins+1);
	}

	else
		log_printf(ERROR, "Bin type %d is not yet implemented for Tally %s",
															_tally_name, type);

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
				 "the bins have not yet been created", _tally_name);

	/* Allocate memory for the bin centers array */
	_centers = new double[_num_bins];

	/* Loop over all bins and find the midpoint between edges */
	for (int i=0; i < _num_bins; i++)
		_centers[i] = (_edges[i] + _edges[i+1]) / 2.0;

	return;
}


/**
 * Tallies a weight for a sample
 * @param neutron the neutron we are tallying
 * @param weight the weight to increment tally by
 */
void Tally::tally(neutron* neutron, double weight) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Tally %s since "
				 "the bins have not yet been created", _tally_name);
	if (_num_batches == 0)
		 log_printf(ERROR, "Cannot tally samples in Tally %s since "
				 "batches have not yet been created", _tally_name);

	int bin_index = getBinIndex(neutron->_old_energy);

	if (bin_index >= 0 && bin_index < _num_bins) 
		_tallies[neutron->_batch_num][bin_index] += weight;

	return;
}


/**
 * Computes average, variance, standard deviation and relative error for each
 * bin over the set of batches
 */
void Tally::computeBatchStatistics() {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot compute batch statistics for Tally %s"
			" since the number of batches has not been set yet", _tally_name);

	double s1 = 0.0;
	double s2 = 0.0;

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {

		s1 = 0.0;
		s2 = 0.0;

		/* Initialize statistics to zero */
		_batch_mu[i] = 0.0;
		_batch_variance[i] = 0.0;
		_batch_std_dev[i] = 0.0;
		_batch_rel_err[i] = 0.0;

		/* Accumulate in s1, s2 counters */
		for (int j=0; j < _num_batches; j++) {
			s1 += _tallies[j][i];
			s2 += _tallies[j][i] * _tallies[j][i];
		}

		/* Compute batch average */
		_batch_mu[i] = s1 / _num_batches;

		/* Compute batch variance */
		_batch_variance[i] = (1.0 / (double(_num_batches) - 1.0)) *
				(s2 / double(_num_batches) - (_batch_mu[i]*_batch_mu[i]));

		_batch_std_dev[i] = sqrt(_batch_variance[i]);
		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
	}

	_computed_statistics = true;

	return;
}


/**
 * Computes average, variance, standard deviation and relative error for each
 * bin over the set of batches. This method first scales each bin value by
 * a scaling factor
 * @param scale_factor the factor to scale each bin value by
 */
void Tally::computeScaledBatchStatistics(double scale_factor) {

	if (_num_batches == 0)
		 log_printf(ERROR, "Cannot compute batch statistics for Tally %s since"
				 " then number of batches has not yet been set", _tally_name);

	double s1 = 0.0;
	double s2 = 0.0;

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {

		s1 = 0.0;
		s2 = 0.0;

		/* Initialize statistics to zero */
		_batch_mu[i] = 0.0;
		_batch_variance[i] = 0.0;
		_batch_std_dev[i] = 0.0;
		_batch_rel_err[i] = 0.0;

		/* Accumulate in s1, s2 counters */
		for (int j=0; j < _num_batches; j++) {
			s1 += _tallies[j][i] / scale_factor;
			s2 += (_tallies[j][i] / scale_factor) *
						(_tallies[j][i] / scale_factor);
		}

		/* Compute batch average */
		_batch_mu[i] = s1 / _num_batches;

		/* Compute batch variance */
		_batch_variance[i] = (1.0 / (double(_num_batches) - 1.0)) *
				(s2 / double(_num_batches) - (_batch_mu[i]*_batch_mu[i]));

		_batch_std_dev[i] = sqrt(_batch_variance[i]);
		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
	}

	_computed_statistics = true;

	return;
}


/**
 * Outputs the batch statistics (if they have been computed) to an
 * ASCII file
 * @param filename the output filename
 */
void Tally::outputBatchStatistics(const char* filename) {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot output batch statistics for Tally %s "
				"since the batches have not yet been generated", _tally_name);

	if (!_computed_statistics)
		log_printf(ERROR, "Cannot output batch statistics for Tally %s "
				"since statistics have not yet been computed", _tally_name);

	/* Create output file */
	FILE* output_file;
	output_file = fopen(filename, "w");

	/* Print header to output file */
	fprintf(output_file, "Batch-based tally statistics for PINSPEC\n");

	if (_tally_type == COLLISION_RATE)
		fprintf(output_file, "Tally type: COLLISION_RATE Rate\n");
	else if (_tally_type == FLUX)
		fprintf(output_file, "Tally type: Flux\n");
	else if (_tally_type == ELASTIC_RATE)
		fprintf(output_file, "Tally type: ELASTIC_RATE Scattering Reaction"
																" Rate\n");
	else if (_tally_type == ABSORPTION_RATE)
		fprintf(output_file, "Tally type: ABSORPTION_RATE Reaction Rate\n");
	else if (_tally_type == CAPTURE_RATE)
		fprintf(output_file, "Tally type: CAPTURE_RATE Reaction Rate\n");
	else if (_tally_type == FISSION_RATE)
		fprintf(output_file, "Tally type: FISSION_RATE Reaction Rate\n");
	else if (_tally_type == TRANSPORT_RATE)
		fprintf(output_file, "Tally type: TRANSPORT_RATE Reaction Rate\n");
	else if (_tally_type == DIFFUSION_RATE)
		fprintf(output_file, "Tally type: DIFFUSION_RATE Reaction Rate\n");
	else if (_tally_type == LEAKAGE_RATE)
		fprintf(output_file, "Tally type: LEAKAGE_RATE Rate\n");

	if (_tally_domain == ISOTOPE)
		fprintf(output_file, "Tally Domain: Isotope\n");
	else if (_tally_domain == MATERIAL)
		fprintf(output_file, "Tally Domain: Material\n");
	else if (_tally_domain == REGION)
		fprintf(output_file, "Tally Domain: Region\n");
	else
		fprintf(output_file, "Tally Domain: Geometry\n");


	if (_bin_spacing == EQUAL)
		fprintf(output_file, "Equally spaced bins with width = %d\n", 
														_bin_spacing);
	else if (_bin_spacing == LOGARITHMIC)
		fprintf(output_file, "Logarithmically spaced bins\n");
	else if (_bin_spacing == OTHER)
		fprintf(output_file, "User-defined bins\n");

	fprintf(output_file, "# batches: %d\t, # bins: %d\n", 
									_num_batches, _num_bins);
	fprintf(output_file, "Bin center, Mu, Variance, Std Dev, Rel Err\n");

	/* Loop over each bin and print mu, var, std dev and rel err */
	for (int i=0; i < _num_bins; i++) {
		fprintf(output_file, "%1.10f, %1.10f, %1.10f, %1.10f, %1.10f\n",
				_centers[i], _batch_mu[i], _batch_variance[i], 
						_batch_std_dev[i], _batch_rel_err[i]);
	}

	fclose(output_file);

	return;
}


/******************************************************************************/
/*********************** Collision Rate Tallying Methods **********************/
/******************************************************************************/

void IsotopeCollisionRateTally::tally(neutron* neutron) {
	Tally::tally(neutron, 1.0);
	return;
}


void MaterialCollisionRateTally::tally(neutron* neutron) {
	Tally::tally(neutron, 1.0);
	return;
}


void RegionCollisionRateTally::tally(neutron* neutron) {
	Tally::tally(neutron, 1.0);
	return;
}


void GeometryCollisionRateTally::tally(neutron* neutron) {
	Tally::tally(neutron, 1.0);
	return;
}

/******************************************************************************/
/************************ Elastic Rate Tallying Methods ***********************/
/******************************************************************************/
//FIXME

void IsotopeElasticRateTally::tally(neutron* neutron) {
	double weight = _isotope->getElasticXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void MaterialElasticRateTally::tally(neutron* neutron) {
	double weight = _material->getElasticMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void RegionElasticRateTally::tally(neutron* neutron) {
	double weight = _region->getElasticMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void GeometryElasticRateTally::tally(neutron* neutron) {
	double weight = neutron->_region->getElasticMicroXS(neutron->_old_energy)
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}

/******************************************************************************/
/********************** Absorption Rate Tallying Methods **********************/
/******************************************************************************/

void IsotopeAbsorptionRateTally::tally(neutron* neutron) {
	double weight = _isotope->getAbsorptionXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
//	log_printf(NORMAL, "Tallied %f for neutron old energy %f with total xs %f for isotope %s absorption rate", weight, neutron->_old_energy, neutron->_total_xs, _tally_name);
	return;
}


void MaterialAbsorptionRateTally::tally(neutron* neutron) {
	double weight = _material->getAbsorptionMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void RegionAbsorptionRateTally::tally(neutron* neutron) {
	double weight = _region->getAbsorptionMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void GeometryAbsorptionRateTally::tally(neutron* neutron) {
	double weight = neutron->_region->getAbsorptionMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}

/******************************************************************************/
/************************ Capture Rate Tallying Methods ***********************/
/******************************************************************************/

void IsotopeCaptureRateTally::tally(neutron* neutron) {
	double weight = _isotope->getCaptureXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void MaterialCaptureRateTally::tally(neutron* neutron) {
	double weight = _material->getCaptureMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void RegionCaptureRateTally::tally(neutron* neutron) {
	double weight = _region->getCaptureMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void GeometryCaptureRateTally::tally(neutron* neutron) {
	double weight = neutron->_region->getCaptureMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


/******************************************************************************/
/*********************** Fission Rate Tallying Methods ************************/
/******************************************************************************/

void IsotopeFissionRateTally::tally(neutron* neutron) {
	double weight = _isotope->getFissionXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void MaterialFissionRateTally::tally(neutron* neutron) {
	double weight = _material->getFissionMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void RegionFissionRateTally::tally(neutron* neutron) {
	double weight = _region->getFissionMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}

void GeometryFissionRateTally::tally(neutron* neutron) {
	double weight = neutron->_region->getFissionMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


/******************************************************************************/
/********************** Transport Rate Tallying Methods ***********************/
/******************************************************************************/

void IsotopeTransportRateTally::tally(neutron* neutron) {
	double weight = _isotope->getTransportXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void MaterialTransportRateTally::tally(neutron* neutron) {
	double weight = _material->getTransportMicroXS(neutron->_old_energy) 
														/neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void RegionTransportRateTally::tally(neutron* neutron) {
	double weight = _region->getTransportMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void GeometryTransportRateTally::tally(neutron* neutron) {
	double weight = neutron->_region->getTransportMicroXS(neutron->_old_energy) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}

/******************************************************************************/
/********************** Diffusion Rate Tallying Methods ***********************/
/******************************************************************************/

void IsotopeDiffusionRateTally::tally(neutron* neutron) {
	double weight = 1.0 / (3.0 * _isotope->getTransportXS(neutron->_old_energy))  
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void MaterialDiffusionRateTally::tally(neutron* neutron) {
	double weight = 1.0 / (3.0 * _material->getTransportMicroXS(neutron->_old_energy)) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void RegionDiffusionRateTally::tally(neutron* neutron) {
	double weight = 1.0 / (3.0 * _region->getTransportMicroXS(neutron->_old_energy)) 
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


void GeometryDiffusionRateTally::tally(neutron* neutron) {
	double weight = 1.0 / (3.0 * neutron->_region->getTransportMicroXS(neutron->_old_energy))
														/ neutron->_total_xs;
	Tally::tally(neutron, weight);
	return;
}


/******************************************************************************/
/*********************** Leakage Rate Tallying Methods ************************/
/******************************************************************************/
//FIXME

void MaterialLeakageRateTally::tally(neutron* neutron) {
	double weight = _material->getBucklingSquared() / 
					(3.0 * _material->getTransportMicroXS(neutron->_old_energy) * 
												neutron->_total_xs);
	Tally::tally(neutron, weight);
	return;
}


void RegionLeakageRateTally::tally(neutron* neutron) {
	double weight = _region->getBucklingSquared() / 
					(3.0 * _region->getTransportMicroXS(neutron->_old_energy) * 
												neutron->_total_xs);
	Tally::tally(neutron, weight);
	return;
}

void GeometryLeakageRateTally::tally(neutron* neutron) {
	double weight = neutron->_region->getBucklingSquared() / 
					(3.0 * neutron->_region->getTransportMicroXS(neutron->_old_energy) * 
												neutron->_total_xs);
	Tally::tally(neutron, weight);
	return;
}


/******************************************************************************/
/************************** Flux Tallying Methods *****************************/
/******************************************************************************/

void MaterialFluxTally::tally(neutron* neutron) {
	Tally::tally(neutron, double(1.0 / neutron->_total_xs));
	return;
}


void RegionFluxTally::tally(neutron* neutron) {
	Tally::tally(neutron, double(1.0 / neutron->_total_xs));
	return;
}


void GeometryFluxTally::tally(neutron* neutron) {
	Tally::tally(neutron, double(1.0 / neutron->_total_xs));
	return;
}
