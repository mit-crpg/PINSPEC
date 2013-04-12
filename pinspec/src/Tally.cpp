#include "Tally.h"


/**
 * @brief Tally constructor.
 * @details Assigns a default number of batches (0) and tally bins (0).
 */
Tally::Tally(const char* tally_name) {

    _tally_name = (char*)tally_name;
    _trigger_type = NONE;
    _trigger_precision = std::numeric_limits<float>::max();

    /* Sets the default delta between bins to zero */
    _bin_delta = 0;

    /* Sets the default for batch statistics */
    _num_batches = 0;
    _num_bins = 0;
    _computed_statistics = false;
}


/**
 * @brief Tally destructor deletes memory for tallies, number of tallies,
 *        bin centers and bin edges (if they have been created).
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
 * @brief Returns the name of the tally.
 * @return the name of the tally
 * @return the Tally's name
 */
char* Tally::getTallyName() {
    return _tally_name;
}


/**
 * @brief Returns the number of tally bins.
 * @return the number of bins
 */
int Tally::getNumBins() {
    return _num_bins;
}


/**
 * @brief Returns a double array of bin edge values.
 * @return array of bin edge values
 */
double* Tally::getBinEdges() {
    if (_num_bins == 0)
        log_printf(ERROR, "Cannot return bin edges for Tally %s since "
		   "the bins have not yet been created", _tally_name);

    return _edges;
}


/**
 * @brief Returns a double array of bin center values.
 * @return array of bin center values
 */
double* Tally::getBinCenters() {
    if (_num_bins == 0)
        log_printf(ERROR, "Cannot return bin centers for Tally %s since the "
		 "centers have not yet been created", _tally_name);

    return _centers;
}


/**
 * @brief Returns the delta spacing between bins. NOTE: this value is only 
 *        non-zero for EQUAL and LOGARITHMIC bin types.
 * @return the spacing between bins
 */
double Tally::getBinDelta() {
    return _bin_delta;
}


/**
 * @brief Returns the delta spacing between the bin edges sandwiching
 *        a given sample value.
 * @param sample the value of interest
 * @return the spacing between bin edges
 */
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
 * @brief Returns the bin spacing type (EQUAL, LOGARITHMIC, OTHER).
 * @return the bin spacing type
 */
binSpacingType Tally::getBinSpacingType() {
    return _bin_spacing;
}


/**
 * @brief Returns the type of tally for these bins (ISOTOPE, MATERIAL, REGION).
 * @return the tally type
 */
tallyDomainType Tally::getTallyDomainType() {
    return _tally_domain;
}


/**
 * @brief Returns the type of tally for these bins (FLUX, COLLISION_RATE, etc).
 * @return the tally type
 */
tallyType Tally::getTallyType() {
    return _tally_type;
}


/**
 * @brief Returns a double array of the tallies within each bin.
 * @return an array of
 */
double** Tally::getTallies() {
    if (_num_bins == 0)
        log_printf(ERROR, "Cannot return tallies for Tally %s since the "
		   "bins have not yet been created", _tally_name);

    return _tallies;
}


/**
 * @brief Returns a specific tally for a specific bin and batch.
 * @param batch_num the batch of interest
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
 * @brief Returns the maximum tally value among all bins and batches.
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
 * @brief Returns the maximum tally value among all bins and batches.
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
 * @brief Returns the number of batches for this tally.
 * @return the number of batches
 */
int Tally::getNumBatches() {
    return _num_batches;
}


/**
 * @brief Returns a pointer to an array of tally batch averages if they 
 *        have been computed.
 * @return a double array of batch averages for each bin
 */
double* Tally::getBatchMu() {

    if (!_computed_statistics)
        log_printf(ERROR, "Statistics have not yet been computed for "
		   "Tally %s so batch mu cannot be returned", _tally_name);

    return _batch_mu;
}


/**
 * @brief Returns a pointer to an array of tally batch variances if they 
 *        have been computed.
 * @return a double array of batch variances for each bin
 */
double* Tally::getBatchVariance() {

    if (!_computed_statistics)
        log_printf(ERROR, "Statistics have not yet been computed for "
		   "Tally %s so batch variance cannot be returned", _tally_name);

    return _batch_variance;
}


/**
 * @brief Returns a pointer to an array of tally batch standard deviations if 
 *        they have been computed.
 * @return a double array of batch standard deviations for each bin
 */
double* Tally::getBatchStdDev() {

    if (!_computed_statistics)
        log_printf(ERROR, "Statistics have not yet been computed for "
		   "Tally %s so batch std dev cannot be returned", _tally_name);

    return _batch_std_dev;
}


/**
 * @brief Returns a pointer to an array of tally batch relative errors if they 
 *        have been computed.
 * @return a double array of batch relative errors for each bin
 */
double* Tally::getBatchRelativeError() {

    if (!_computed_statistics)
        log_printf(ERROR, "Statistics have not yet been computed for "
		   "Tally %s so batch relative error cannot be returned", 
		   _tally_name);

    return _batch_rel_err;
}


/**
 * @brief Returns the maximum average tally over batches.
 * @return the maximum average tally
 */
double Tally::getMaxMu() {

    double max_mu = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_mu[i] > max_mu)
            max_mu = _batch_mu[i];
    }

    return max_mu;
}


/**
 * @brief Returns the maximum tally variance over batches.
 * @return the maximum tally variance
 */
double Tally::getMaxVariance() {

    double max_variance = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_variance[i] > max_variance)
            max_variance = _batch_variance[i];
    }

    return max_variance;
}


/**
 * @brief Returns the maximum tally relative error.
 * @return the maximum relative error
 */
double Tally::getMaxRelErr() {

    double max_rel_err = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_rel_err[i] > max_rel_err)
            max_rel_err = _batch_rel_err[i];
    }

    return max_rel_err;
}


/**
 * @brief Returns the maximum tally standard deviatoin.
 * @return the maximum standard deviation
 */
double Tally::getMaxStdDev() {

    double max_std_dev = 0.0;

    for (int i=0; i < _num_bins; i++) {
        if (_batch_std_dev[i] > max_std_dev)
            max_std_dev = _batch_std_dev[i];
    }

    return max_std_dev;
}


/**
 * @brief Returns the trigger precision for this tally.
 * @details The trigger precision is the threshold value which must be met 
 *          for all of this tally's values before the PINSPEC simulation will 
 *          complete.
 * @return the trigger precision
 */
float Tally::getTriggerPrecision() {
    return _trigger_precision;
}


/**
 * @brief Returns the precision trigger type (VARIANCE, STANDARD_DEVIATION,
 *        RELATIVE_ERROR, or NONE).
 * @return the trigger precision type.
 */
triggerType Tally::getTriggerType() {
    return _trigger_type;
}


/**
 * @brief Returns whether or not the tally has computed batch statistics.
 * @return true if the tally has computed batch statistics; otherwise false
 */
bool Tally::hasComputedBatchStatistics() {    
    return _computed_statistics;
}


/**
 * @brief Returns whether or not the tally precision meets the 
 *        precision trigger threshold, if a trigger exists.
 * @return true if the precision meets the threshold; otherwise false
 */
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



/**
 * @brief This method fills an array with the tally bin edges. 
 * @details This method is a utility function for users to access PINSPEC data 
 *          in Python. The method prototype may seem to require two arguments 
 *          - the array to fill and the number of tally bins - but in Python 
 *          the user must only supply the data array as follows:
 *
 * @code
 *        num_bins = tally.getNumBins()
 *        bin_edges = numpy.zeros(num_bins)
 *        tally.retrieveTallyEdges(bin_edges)
 * @endcode
 *
 * @param data the data array to fill with bin edge values
 * @param num_bins the number of tally bins
 */
void Tally::retrieveTallyEdges(double* data, int num_bins) {

    /* Load all tally bin centers into array */
    for (int i=0; i < _num_bins+1; i++)
        data[i] = _edges[i];
}


/**
 * @brief This method fills an array with the tally bin centers. 
 * @details This method is a utility function for users to access PINSPEC data 
 *          in Python. The method prototype may seem to require two arguments 
 *          - the array to fill and the number of tally bins - but in Python 
 *          the user must only supply the data array as follows:
 *
 * @code
 *        num_bins = tally.getNumBins()
 *        bin_centers = numpy.zeros(num_bins)
 *        tally.retrieveTallyCenters(bin_centers)
 * @endcode
 *
 * @param data the data array to fill with bin center values
 * @param num_bins the number of tally bins
 */
void Tally::retrieveTallyCenters(double* data, int num_bins){

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve bin centers for Tally %s since"
                      " it has not yet computed batch statistics", _tally_name);

    /* Load all tally bin centers into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _centers[i];
}


/**
 * @brief This method fills an array with the average tally values. 
 * @details This method is a utility function for users to access PINSPEC data 
 *          in Python. The method prototype may seem to require two arguments 
 *          - the array to fill and the number of tally bins - but in Python 
 *          the user must only supply the data array as follows:
 *
 * @code
 *        num_bins = tally.getNumBins()
 *        tally_averages = numpy.zeros(num_bins)
 *        tally.retrieveTallyMu(tally_averages)
 * @endcode
 *
 * @param data the data array to fill with tally average values
 * @param num_bins the number of tally bins
 */
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


/**
 * @brief This method fills an array with the tally variances. 
 * @details This method is a utility function for users to access PINSPEC data 
 *          in Python. The method prototype may seem to require two arguments 
 *          - the array to fill and the number of tally bins - but in Python 
 *          the user must only supply the data array as follows:
 *
 * @code
 *        num_bins = tally.getNumBins()
 *        variances = numpy.zeros(num_bins)
 *        tally.retrieveTallyVariance(variances)
 * @endcode
 *
 * @param data the data array to fill with tally variances
 * @param num_bins the number of tally bins
 */
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


/**
 * @brief This method fills an array with the tally standard deviations. 
 * @details This method is a utility function for users to access PINSPEC data 
 *          in Python. The method prototype may seem to require two arguments - 
 *          the array to fill and the number of tally bins - but in Python the 
 *          user must only supply the data array as follows:
 *
 * @code
 *        num_bins = tally.getNumBins()
 *        std_dev = numpy.zeros(num_bins)
 *        tally.retrieveTallyStdDev(std_dev)
 * @endcode
 *
 * @param data the data array to fill with tally standard deviations
 * @param num_bins the number of tally bins
 */
void Tally::retrieveTallyStdDev(double* data, int num_bins) {

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve tally std. dev. for Tally %s "
		   "since it has not yet computed batch statistics", 
		   _tally_name);

    if (_num_batches == 0)
        log_printf(ERROR, "Unable to retrieve tally std. dev. for Tally %s "
		   "since it does not know how many batches it should tally", 
		   _tally_name);

    /* Load all tally standard deviations into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _batch_std_dev[i];
}


/**
 * @brief This method fills an array with the tally relative errors. 
 * @details This method is a utility function for users to access PINSPEC data 
 *          in Python. The method prototype may seem to require two arguments - 
 *          the array to fill and the number of tally bins - but in Python the 
 *          user must only supply the data array as follows:
 *
 * @code
 *        num_bins = tally.getNumBins()
 *        rel_err = numpy.zeros(num_bins)
 *        tally.retrieveTallyRelErr(rel_err)
 * @endcode
 *
 * @param data the data array to fill with tally relative errors
 * @param num_bins the number of tally bins
 */
void Tally::retrieveTallyRelErr(double* data, int num_bins) {

    if (!_computed_statistics)
        log_printf(ERROR, "Unable to retrieve tally rel. err. for Tally %s "
		   "since it has not yet computed batch statistics", 
		   _tally_name);

    if (_num_batches == 0)
        log_printf(ERROR, "Unable to retrieve tally rel. err. for Tally %s "	
		   "since it does not know how many batches it should tally", 
		   _tally_name);

    /* Load all tally variances into array */
    for (int i=0; i < _num_bins; i++)
        data[i] = _batch_rel_err[i];
}


/**
 * @brief Set the bin spacing type for this Tally (EQUAL, LOGARITHMIC, OTHER).
 * @param type the bin spacing type
 */
void Tally::setBinSpacingType(binSpacingType type) {
    _bin_spacing = type;
}


/**
 * @brief Set the tally domain type (MATERIAL, REGION, etc.).
 * @param type the tally domain type
 */
void Tally::setTallyDomainType(tallyDomainType type) {
    _tally_domain = type;
}


/**
 * @brief Set the tally type (FLUX, CAPTURE_RATE, etc.).
 * @param type the tally type
 */
void Tally::setTallyType(tallyType type) {
    _tally_type = type;
}



/**
 * @brief Set a user-defined double array of bin edge values.
 * @details This method is intended to allow PINSPEC users to set tally bin 
 *          edges through Python. Although this method prototype seems to 
 *          require two arguments - an array of bin edges and the number of 
 *          edges - in Python the user only needs to provide the edges array 
 *          as follows:
 *
 * @code
 *          edges = numpy.array([0.1, 1., 5., 25., 100.])
 *          tally.setBinEdges(edges)
 * @endcode
 *
 * @param edges the array of bin edges
 * @param num_edges the number of bin edges
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


/**
 * @brief Sets a precision trigger for this tally.
 * @details By setting a precision trigger, the user instructs a PINSEPC
 *          simulation to continue running until all tallies meet the
 *          precision trigger threshold.
 * @param trigger_type the precision trigger type (VARIANCE, etc)
 * @param precision the threshold for the precision trigger
 */
void Tally::setPrecisionTrigger(triggerType trigger_type, float precision) {

    if (precision < 0.0)
	log_printf(ERROR, "Unable to set a negative trigger precision of %f"
		   " for tally %s", precision, _tally_name);

    _trigger_type = trigger_type;
    _trigger_precision = precision;
    return;
}


/**
 * @brief Set the number of batches for this Tally. 
 * @details This method also allocates memory for the tallies and batch 
 *          statistics arrays.
 * @param num_batches the number of batches
 */
void Tally::setNumBatches(int num_batches) {

    /* Clean up memory from old arrays of batch statistics */
    if (_num_batches != 0) {
      //        delete [] _tallies;
      //	delete [] _batch_mu;
      //	delete [] _batch_variance;
      //	delete [] _batch_std_dev;
      //	delete [] _batch_rel_err;
    }

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


/**
 * @brief Increments the number of batches for this tally.
 * @param num_batches number of batches to add to the total number of batches.
 */
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
 * @brief Generate edges between bins defined by a start and end point.
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

    /* Generate arrays for a default of 1 set of batch statistics */
    setNumBatches(1);

    return;
}


/**
 * @brief Compute the center points between bin edges for this Tally's bins.
 * @details This method populates a private class attribute array with the
 *          bin center values.
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
 * @brief This method tallies a particular weight for a neutron.
 * @details The method determines which tally bin to use based on the
 *          neutron's energy.
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

    if (weight < 0.0)
        log_printf(NORMAL, "weight = %f", weight);

    if (bin_index >= 0 && bin_index < _num_bins) 
        _tallies[neutron->_batch_num][bin_index] += weight;

    return;
}


/**
 * @brief Computes average, variance, standard deviation and relative error 
 *        for each bin over the set of batches.
 * @details This method populates private class attribute arrays with the
 *          batch statistics.
 */
void Tally::computeBatchStatistics() {

    if (_num_batches == 0)
        log_printf(ERROR, "Cannot compute batch statistics for Tally %s"
		   " since the number of batches has not been set yet", 
		   _tally_name);

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
 * @brief Computes average, variance, standard deviation and relative error 
 *        for each bin over the set of batches. 
 * @details This method scales each bin value by a scaling factor.
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
				(s2 / double(_num_batches) - 
				 (_batch_mu[i]*_batch_mu[i]));

	_batch_std_dev[i] = sqrt(_batch_variance[i]);
	_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
    }
    
    _computed_statistics = true;
    
    return;
}


/**
 * @brief Divide each tally by the maximum tally value.
 */
void Tally::normalizeBatchMu() {

    if (_num_bins == 0)
        log_printf(ERROR, "Cannot normalize batch mu for Tally %s since its "
		   " bins have not yet been created", _tally_name);
    if (!_computed_statistics)
      log_printf(ERROR, "Cannot normalize batch mu for Tally %s since it "
		 " has not yet computed batch statistics", _tally_name);

    double max_mu = getMaxMu();

    /* Divide each tally by maximum tally value */
    for (int n=0; n < _num_bins; n++)
        _batch_mu[n] /= max_mu;

    return;
}

/**
 * @brief Outputs the batch statistics (if they have been computed) to an
 *        ASCII file.
 * @param filename the output filename (optional)
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
    fprintf(output_file, "Tally name: %s\n", _tally_name);

    if (_tally_type == COLLISION_RATE)
        fprintf(output_file, "Tally type: COLLISION_RATE Rate\n");
    else if (_tally_type == FLUX)
        fprintf(output_file, "Tally type: Flux\n");
    else if (_tally_type == ELASTIC_RATE)
        fprintf(output_file, "Tally type: ELASTIC_RATE Scattering Reaction "
		"Rate\n");
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
        fprintf(output_file, "Equally spaced bins with width = %f\n", 
		_bin_delta);
    else if (_bin_spacing == LOGARITHMIC)
        fprintf(output_file, "Logarithmically spaced bins with"
		" width = %f\n", _bin_delta);
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


/**
 * @brief Print the tally values to the screen. 
 * @details This method will print the tally results, including bin edges
 *          tally averages, and if a user requests, even tally batch
 *          statistics.
 * @param uncertainties a boolean representing whether or not to print 
 *        tally statistics
 */
void Tally::printTallies(bool uncertainties) {

    log_printf(HEADER, "Batch Statistics for Tally %s", _tally_name);

    /* Create a tally statistics table header */
    std::stringstream title;
    title << std::string(7, ' ') << "Energy Band" << std::string(9, ' ');
    title << "   Mu   ";

    if (uncertainties) {
        title << std::string(2, ' ') << "Variance";
        title << std::string(2, ' ') << "Std. Dev.";
        title << std::string(1, ' ') << "Rel. Err.";
    }

    log_printf(SEPARATOR, "");
    log_printf(RESULT, title.str().c_str());
    log_printf(SEPARATOR, "");

    char lower_bound[16];
    char upper_bound[16];
    char mu[12];
    char variance[16];
    char std_dev[16];
    char rel_err [16];

    /* Loop over each reaction rate and construct a string with the reate
     * and the user-specified statistics to print to the shell */
    for (int i=0; i < _num_bins; i++) {

        std::stringstream entry;

        /* Format lower bound of bin energy interval */
        if (_edges[i] == 0.0)
            sprintf(lower_bound, "%7.2f", _edges[i]);
        else if (_edges[i] > 0.0 && _edges[i] < 1E-2)
            sprintf(lower_bound, "%7.1E", _edges[i]);
        else if (_edges[i] >= 1E-2 && _edges[i] < 1E4)
            sprintf(lower_bound, "%7.2f", _edges[i]);
        else
            sprintf(lower_bound, "%7.1E", _edges[i]);

        /* Format upper bound of bin energy interval */
        if (_edges[i+1] == 0.0)
            sprintf(upper_bound, "%7.2f", _edges[i+1]);
        else if (_edges[i+1] > 0.0 && _edges[i+1] < 1E-2)
            sprintf(upper_bound, "%7.1E", _edges[i+1]);
        else if (_edges[i+1] > 1E-2 && _edges[i+1] < 1E4)
            sprintf(upper_bound, "%7.2f", _edges[i+1]);
        else
            sprintf(upper_bound, "%7.1E", _edges[i+1]);

        /* Format batch average */
        if (_batch_mu[i] < 1E-2)
            sprintf(mu, "%8.2E", _batch_mu[i]);
        else if (_batch_mu[i] >= 1E-2 && _batch_mu[i] < 10.)
            sprintf(mu, "%8.6f", _batch_mu[i]);
        else if (_batch_mu[i] >= 10. &&  _batch_mu[i] < 1E2)
            sprintf(mu, "%8.5f", _batch_mu[i]);
        else if (_batch_mu[i] >= 1E2 && _batch_mu[i] < 1E3)
            sprintf(mu, "%8.4f", _batch_mu[i]);
        else if (_batch_mu[i] >= 1E3 && _batch_mu[i] < 1E4)
            sprintf(mu, "%8.3f", _batch_mu[i]);
        else if (_batch_mu[i] >= 1E4 && _batch_mu[i] < 1E5)
            sprintf(mu, "%8.2f", _batch_mu[i]);
        else if (_batch_mu[i] >= 1E5 && _batch_mu[i] < 1E6)
            sprintf(mu, "%8.2f", _batch_mu[i]);
        else
            sprintf(mu, "%8.2E", _batch_mu[i]);


        entry << "[ " << lower_bound << " - " << upper_bound << " eV ]:  ";
        entry << mu;
        
        /* No need to format uncertainties since we can assume they are small
         * and use scientific notation */
        if (uncertainties) {
            sprintf(variance, "%8.2E", _batch_variance[i]);
            sprintf(std_dev, "%8.2E", _batch_std_dev[i]);
            sprintf(rel_err, "%8.2E", _batch_rel_err[i]);
            entry << std::string(2, ' ') << variance;
            entry << std::string(2, ' ') << std_dev;
            entry << std::string(2, ' ') << rel_err;
        }

        log_printf(RESULT, entry.str().c_str());
    }

    log_printf(SEPARATOR, "");
}


/**
 * @brief Creates a new version of this tally with identical data.
 * @details The clone method makes a deep copy of all of this tally's
 *          data and loads it into a new tally class object.
 */
Tally* Tally::clone() {
    
    /* Deep copy all of the parameters attributes to this Tally */
    DerivedTally* tally = new DerivedTally(_tally_name);
    tally->setTallyDomainType(_tally_domain);
    tally->setTallyType(_tally_type);
    tally->setBinSpacingType(_bin_spacing);
    tally->setPrecisionTrigger(_trigger_type, _trigger_precision);

    /* If we have LOGARITHMIC or EQUAL bins, generate them for the new clone */
    if (_bin_spacing == LOGARITHMIC || _bin_spacing == EQUAL)
        tally->generateBinEdges(_edges[0], _edges[_num_bins], 
                                                    _num_bins, _bin_spacing);
    else
        tally->setBinEdges(_edges, _num_bins+1);

    tally->setNumBatches(_num_batches);
    tally->setComputedBatchStatistics(_computed_statistics);

    if (_computed_statistics) {
        double* batch_mu = new double[_num_bins];
        double* batch_variance = new double[_num_bins];
        double* batch_std_dev = new double[_num_bins];
        double* batch_rel_err = new double[_num_bins];

        memcpy(batch_mu, _batch_mu, _num_bins*sizeof(double));
        memcpy(batch_variance, _batch_variance, _num_bins*sizeof(double));
        memcpy(batch_std_dev, _batch_std_dev, _num_bins*sizeof(double));
        memcpy(batch_rel_err, _batch_rel_err, _num_bins*sizeof(double));

        tally->setBatchMu(batch_mu);
        tally->setBatchVariance(batch_variance);
        tally->setBatchStdDev(batch_std_dev);
        tally->setBatchRelErr(batch_rel_err);
    }

    double** tallies = (double**) malloc(sizeof(double*) * _num_batches);
    for (int i=0; i < _num_batches; i++) {
	    tallies[i] = new double[_num_bins];
        memcpy(tallies[i], _tallies[i], _num_bins);
    }
    
    tally->setTallies(tallies);

    return tally;
}


/**
 * @brief Tally addition operator.
 * @details This overloaded addition operator allows the user to
 *          add two tallies with each other, if they have the same number 
 *          of tallies. The creates a new DERIVED type tally, loads it with 
 *          the sum of the tally bin averages for the two tally operands, 
 *          and computes its batch statistics appropriately. This method is 
 *          intended to allow for simple tally arithmetic in Python.
 * @code 
 *          new_tally = tally1 + tally2
 * @endcode
 * @param tally the right operand in the tally summation
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::operator+(Tally* tally) {

    /* If the tallies have not yet computed batch statistics, reject */
    if (!_computed_statistics)
        log_printf(ERROR, "Unable to add tally %s which has not yet computed "
                                            "batch statistics", _tally_name);
    else if (!tally->hasComputedBatchStatistics())
        log_printf(ERROR, "Unable to add tally %s which has not yet computed "
                                    "batch statistics", tally->getTallyName());

    /* Check that the tally bin numbers and edges match appropriately */
    else if (_num_bins != 1 && tally->getNumBins() != 1) {

        /* If the two tallies have different numbers of bins, reject */
        if (_num_bins != tally->getNumBins())
            log_printf(ERROR, "Unable to add tally %s with %d bins to tally %s"
                                " with %d bins", _tally_name, _num_bins,
                                tally->getTallyName(), tally->getNumBins());

        /* Check that all bin edges match */
        double* edges = new double[_num_bins+1];
        tally->retrieveTallyEdges(edges, _num_bins);

        /* Loop over all edges */
        for (int i=0; i < _num_bins+1; i++) {
            if (_edges[i] != edges[i])
                log_printf(ERROR, "Unable to add tally %s with bin edge %d to "
                            "tally %s with bin edge %d", 
                            _tally_name, _edges[i],
                           tally->getTallyName(), edges[i]);
        }
    }

    /* Create a new derived type tally and initialize it */        
    const char* tally_name = (std::string(_tally_name) + 
                                    " + " + tally->getTallyName()).c_str(); 
    DerivedTally* new_tally = new DerivedTally(tally_name);

    new_tally->setBinSpacingType(_bin_spacing);
    new_tally->setBinEdges(_edges, _num_bins+1);

    /* Initialize new arrays for the derived tallies statistics */
    int max_num_bins = std::max(_num_bins, tally->getNumBins());
    double* new_mu = new double[max_num_bins];
    double* new_variance = new double[max_num_bins];
    double* new_std_dev = new double[max_num_bins];
    double* new_rel_err = new double[max_num_bins];

    /* Retrieve statistics metrics from RHS tally */
    double* mu = new double[tally->getNumBins()];
    double* variance = new double[tally->getNumBins()];

    tally->retrieveTallyMu(mu, tally->getNumBins());
    tally->retrieveTallyVariance(variance ,tally->getNumBins());

    double *mu1, *mu2;
    double *variance1, *variance2;

    /* Check if the LHS tally (this one) only has one bin */
    if (_num_bins == max_num_bins) {
        mu1 = _batch_mu;
        variance1 = _batch_variance;
    }
    else {
        mu1 = new double[max_num_bins];
        variance1 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu1[i] = _batch_mu[0];
            variance1[i] = _batch_variance[0];
        }
    }

    /* Check if the RHS tally only has one bin */        
    if (tally->getNumBins() == max_num_bins) {
        mu2 = mu;
        variance2 = variance;
    }
    else {
        mu2 = new double[max_num_bins];
        variance2 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu2[i] = mu[0];
            variance2[i] = variance[0];
        }
    }

    /* Compute the derived tallies batch mu */
    for (int i=0; i < max_num_bins; i++) {
        new_mu[i] = mu1[i] + mu2[i];
        new_variance[i] = variance1[i] + variance2[i];
        new_std_dev[i] = sqrt(new_variance[i]);
        new_mu[i] = new_std_dev[i] / new_mu[i];
    }

    new_tally->setNumBatches(1);
    new_tally->setBatchMu(new_mu);
    new_tally->setBatchVariance(new_variance);
    new_tally->setBatchStdDev(new_std_dev);
    new_tally->setBatchRelErr(new_rel_err);
    new_tally->setComputedBatchStatistics(true); 

    return new_tally;
}


/**
 * @brief Tally subtraction operator
 * @details This overloaded subtraction operator allows the user to
 *          subtract two tallies with each other, if they have the same number 
 *          of tallies. The creates a new DERIVED type tally, loads it with 
 *          the difference of the tally bin averages for the two tally operands
 *          and computes its batch statistics appropriately. This method is 
 *          intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally1 - tally2
 * @endcode
 *
 * @param tally the right operand in the tally product
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::operator-(Tally* tally) {

    /* If the tallies have not yet computed batch statistics, reject */
    if (!_computed_statistics)
        log_printf(ERROR, "Unable to subtract tally %s which has not yet "
                        "computed batch statistics", _tally_name);
    else if (!tally->hasComputedBatchStatistics())
        log_printf(ERROR, "Unable to subtract tally %s which has not yet "
                        "computed batch statistics", tally->getTallyName());

    /* Check that the tally bin numbers and edges match appropriately */
    else if (_num_bins != 1 && tally->getNumBins() != 1) {

        /* If the two tallies have different numbers of bins, reject */
        if (_num_bins != tally->getNumBins())
            log_printf(ERROR, "Unable to subtract tally %s with %d bins and "
                            " tally %s with %d bins", _tally_name, _num_bins,
                            tally->getTallyName(), tally->getNumBins());

        /* Check that all bin edges match */
        double* edges = new double[_num_bins+1];
        tally->retrieveTallyEdges(edges, _num_bins);

        /* Loop over all edges */
        for (int i=0; i < _num_bins+1; i++) {
            if (_edges[i] != edges[i])
                log_printf(ERROR, "Unable to subtract tally %s with bin edge"
                            " %d and tally %s with bin edge %d", 
                            _tally_name, _edges[i],
                           tally->getTallyName(), edges[i]);
        }
    }

    /* Create a new derived type tally and initialize it */        
    const char* tally_name = (std::string(_tally_name) + 
                                      " - " + tally->getTallyName()).c_str(); 
    DerivedTally* new_tally = new DerivedTally(tally_name);

    new_tally->setBinSpacingType(_bin_spacing);
    new_tally->setBinEdges(_edges, _num_bins+1);

    /* Initialize new arrays for the derived tallies statistics */
    int max_num_bins = std::max(_num_bins, tally->getNumBins());
    double* new_mu = new double[max_num_bins];
    double* new_variance = new double[max_num_bins];
    double* new_std_dev = new double[max_num_bins];
    double* new_rel_err = new double[max_num_bins];

    /* Retrieve statistics metrics from RHS tally */
    double* mu = new double[tally->getNumBins()];
    double* variance = new double[tally->getNumBins()];

    tally->retrieveTallyMu(mu, tally->getNumBins());
    tally->retrieveTallyVariance(variance ,tally->getNumBins());

    double *mu1, *mu2;
    double *variance1, *variance2;

    /* Check if the LHS tally (this one) only has one bin */
    if (_num_bins == max_num_bins) {
        mu1 = _batch_mu;
        variance1 = _batch_variance;
    }
    else {
        mu1 = new double[max_num_bins];
        variance1 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu1[i] = _batch_mu[0];
            variance1[i] = _batch_variance[0];
        }
    }

    /* Check if the RHS tally only has one bin */        
    if (tally->getNumBins() == max_num_bins) {
        mu2 = mu;
        variance2 = variance;
    }
    else {
        mu2 = new double[max_num_bins];
        variance2 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu2[i] = mu[0];
            variance2[i] = variance[0];
        }
    }

    /* Compute the derived tallies batch mu */
    for (int i=0; i < max_num_bins; i++) {
        new_mu[i] = mu1[i] - mu2[i];
        new_variance[i] = variance1[i] + variance2[i];
        new_std_dev[i] = sqrt(new_variance[i]);
        new_mu[i] = new_std_dev[i] / new_mu[i];
    }

    new_tally->setNumBatches(1);
    new_tally->setBatchMu(new_mu);
    new_tally->setBatchVariance(new_variance);
    new_tally->setBatchStdDev(new_std_dev);
    new_tally->setBatchRelErr(new_rel_err);
    new_tally->setComputedBatchStatistics(true); 

    return new_tally;
}


/**
 * @brief Tally multiplication operator.
 * @details This overloaded multiplication operator allows the user to
 *          multiply two tallies with each other, if they have the same number 
 *          of tallies. The creates a new DERIVED type tally, loads it with 
 *          the product of the tally bin averages for the two tally operands, 
 *          and computes its batch statistics appropriately. This method is 
 *          intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally1 * tally2
 * @endcode
 *
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator*(Tally* tally) {

    /* If the tallies have not yet computed batch statistics, reject */
    if (!_computed_statistics)
        log_printf(ERROR, "Unable to multiply tally %s which has not yet "
                        "computed batch statistics", _tally_name);
    else if (!tally->hasComputedBatchStatistics())
        log_printf(ERROR, "Unable to multiply tally %s which has not yet "
                        "computed batch statistics", tally->getTallyName());

    /* Check that the tally bin numbers and edges match appropriately */
    else if (_num_bins != 1 && tally->getNumBins() != 1) {

        /* If the two tallies have different numbers of bins, reject */
        if (_num_bins != tally->getNumBins())
            log_printf(ERROR, "Unable to multiply tally %s with %d bins and "
                            " tally %s with %d bins", _tally_name, _num_bins,
                            tally->getTallyName(), tally->getNumBins());

        /* Check that all bin edges match */
        double* edges = new double[_num_bins+1];
        tally->retrieveTallyEdges(edges, _num_bins);

        /* Loop over all edges */
        for (int i=0; i < _num_bins+1; i++) {
            if (_edges[i] != edges[i])
                log_printf(ERROR, "Unable to multiply tally %s with bin edge"
                            " %d and tally %s with bin edge %d", 
                            _tally_name, _edges[i],
                           tally->getTallyName(), edges[i]);
        }
    }

    /* Create a new derived type tally and initialize it */        
    const char* tally_name = (std::string(_tally_name) + 
                                      " * " + tally->getTallyName()).c_str(); 
    DerivedTally* new_tally = new DerivedTally(tally_name);

    new_tally->setBinSpacingType(_bin_spacing);
    new_tally->setBinEdges(_edges, _num_bins+1);

    /* Initialize new arrays for the derived tallies statistics */
    int max_num_bins = std::max(_num_bins, tally->getNumBins());
    double* new_mu = new double[max_num_bins];
    double* new_variance = new double[max_num_bins];
    double* new_std_dev = new double[max_num_bins];
    double* new_rel_err = new double[max_num_bins];

    /* Retrieve statistics metrics from RHS tally */
    double* mu = new double[tally->getNumBins()];
    double* variance = new double[tally->getNumBins()];

    tally->retrieveTallyMu(mu, tally->getNumBins());
    tally->retrieveTallyVariance(variance ,tally->getNumBins());

    double *mu1, *mu2;
    double *variance1, *variance2;

    /* Check if the LHS tally (this one) only has one bin */
    if (_num_bins == max_num_bins) {
        mu1 = _batch_mu;
        variance1 = _batch_variance;
    }
    else {
        mu1 = new double[max_num_bins];
        variance1 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu1[i] = _batch_mu[0];
            variance1[i] = _batch_variance[0];
        }
    }

    /* Check if the RHS tally only has one bin */        
    if (tally->getNumBins() == max_num_bins) {
        mu2 = mu;
        variance2 = variance;
    }
    else {
        mu2 = new double[max_num_bins];
        variance2 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu2[i] = mu[0];
            variance2[i] = variance[0];
        }
    }

    /* Compute the derived tallies batch mu */
    for (int i=0; i < max_num_bins; i++) {
        new_mu[i] = mu1[i] * mu2[i];
        new_variance[i] = mu1[i]*mu1[i]*variance2[i] + 
                            mu2[i]*mu2[i]*variance1[i] + 
                            variance1[i] * variance2[i];
        new_std_dev[i] = sqrt(new_variance[i]);
        new_mu[i] = new_std_dev[i] / new_mu[i];
    }

    new_tally->setNumBatches(1);
    new_tally->setBatchMu(new_mu);
    new_tally->setBatchVariance(new_variance);
    new_tally->setBatchStdDev(new_std_dev);
    new_tally->setBatchRelErr(new_rel_err);
    new_tally->setComputedBatchStatistics(true); 

    return new_tally;
}


/**
 * @brief Tally division operator.
 * @details This overloaded division operator allows the user to
 *          divide two tallies with each other, if they have the same number 
 *          of tallies. The creates a new DERIVED type tally, loads it with 
 *          the dividend of the tally bin averages for the two tally operands, 
 *          and computes its batch statistics appropriately. This method is 
 *          intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally1 / tally2
 * @endcode
 *
 * @param tally the right operand in the tally division
 * @return a DERIVED type tally with the tally division
 */
DerivedTally* Tally::operator/(Tally* tally) {

    /* If the tallies have not yet computed batch statistics, reject */
    if (!_computed_statistics)
        log_printf(ERROR, "Unable to divide tally %s which has not yet "
                        "computed batch statistics", _tally_name);
    else if (!tally->hasComputedBatchStatistics())
        log_printf(ERROR, "Unable to divide tally %s which has not yet "
                        "computed batch statistics", tally->getTallyName());

    /* Check that the tally bin numbers and edges match appropriately */
    else if (_num_bins != 1 && tally->getNumBins() != 1) {

        /* If the two tallies have different numbers of bins, reject */
        if (_num_bins != tally->getNumBins())
            log_printf(ERROR, "Unable to divide tally %s with %d bins and "
                            " tally %s with %d bins", _tally_name, _num_bins,
                            tally->getTallyName(), tally->getNumBins());

        /* Check that all bin edges match */
        double* edges = new double[_num_bins+1];
        tally->retrieveTallyEdges(edges, _num_bins);

        /* Loop over all edges */
        for (int i=0; i < _num_bins+1; i++) {
            if (_edges[i] != edges[i])
                log_printf(ERROR, "Unable to divide tally %s with bin edge"
                            " %d and tally %s with bin edge %d", 
                            _tally_name, _edges[i],
                           tally->getTallyName(), edges[i]);
        }
    }

    /* Create a new derived type tally and initialize it */        
    const char* tally_name = (std::string(_tally_name) + 
                                       " / " + tally->getTallyName()).c_str(); 
    DerivedTally* new_tally = new DerivedTally(tally_name);

    new_tally->setBinSpacingType(_bin_spacing);
    new_tally->setBinEdges(_edges, _num_bins+1);

    /* Initialize new arrays for the derived tallies statistics */
    int max_num_bins = std::max(_num_bins, tally->getNumBins());
    double* new_mu = new double[max_num_bins];
    double* new_variance = new double[max_num_bins];
    double* new_std_dev = new double[max_num_bins];
    double* new_rel_err = new double[max_num_bins];

    /* Retrieve statistics metrics from RHS tally */
    double* mu = new double[tally->getNumBins()];
    double* variance = new double[tally->getNumBins()];

    tally->retrieveTallyMu(mu, tally->getNumBins());
    tally->retrieveTallyVariance(variance ,tally->getNumBins());

    double *mu1, *mu2;
    double *variance1, *variance2;

    /* Check if the LHS tally (this one) only has one bin */
    if (_num_bins == max_num_bins) {
        mu1 = _batch_mu;
        variance1 = _batch_variance;
    }
    else {
        mu1 = new double[max_num_bins];
        variance1 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu1[i] = _batch_mu[0];
            variance1[i] = _batch_variance[0];
        }
    }

    /* Check if the RHS tally only has one bin */        
    if (tally->getNumBins() == max_num_bins) {
        mu2 = mu;
        variance2 = variance;
    }
    else {
        mu2 = new double[max_num_bins];
        variance2 = new double[max_num_bins];
        for (int i=0; i < max_num_bins; i++) {
            mu2[i] = mu[0];
            variance2[i] = variance[0];
        }
    }

    /* Compute the derived tallies batch mu */
    /* http://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables */        
    for (int i=0; i < max_num_bins; i++) {
        new_mu[i] = (mu1[i] / mu2[i]) 
                    + (mu1[i] / (mu2[i] * mu2[i] * mu2[i])) * variance2[i];
        new_variance[i] = variance1[i] / (mu2[i] * mu2[i]) +
                            ((mu1[i] * mu1[i] * variance2[i]) /
                            (mu2[i] * mu2[i] * mu2[i] * mu2[i]));
        new_std_dev[i] = sqrt(new_variance[i]);
        new_rel_err[i] = new_std_dev[i] / new_mu[i];
    }

    new_tally->setNumBatches(1);
    new_tally->setBatchMu(new_mu);
    new_tally->setBatchVariance(new_variance);
    new_tally->setBatchStdDev(new_std_dev);
    new_tally->setBatchRelErr(new_rel_err);
    new_tally->setComputedBatchStatistics(true); 

    return new_tally;
}


/**
 * @brief Tally addition with a constant operator.
 * @details This overloaded addition operator allows the user to
 *          add a constant to a tally. The creates a new DERIVED type
 *          tally, loads it with the sum of the tally bin averages and 
 *          the constant, and updates its batch statistics appropriately.
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally + int(3)
 * @endcode
 *
 * @param amt the right operand and the constant value to add to the tally
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::operator+(int amt) {
    
    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    for (int i=0; i < _num_bins; i++)
        new_tally->_batch_mu[i] += amt;

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally subtraction with a constant operator.
 * @details This overloaded subtraction operator allows the user to
 *          subtract a constant from a tally. The creates a new DERIVED type
 *          tally, loads it with the difference of the tally bin averages and 
 *          the constant, and updates its batch statistics appropriately.
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally - int(3)
 * @endcode
 *
 * @param amt the right operand and the constant value to subtract from the
 *        tally
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::operator-(const int amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    for (int i=0; i < _num_bins; i++) 
        new_tally->_batch_mu[i] -= amt;

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally multiplication with a constant operator.
 * @details This overloaded subtraction operator allows the user to
 *          multiply a constant with a tally. The creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          constant, and updates its batch statistics appropriately.
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally * int(3)
 * @endcode
 *
 * @param amt the right operand and the constant value to multiply the tally by
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator*(const int amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Multiplying this tally by some amount multiplies the batch mu
     * by a constant and updates the perturbed batch statistics */
    for (int i=0; i < _num_bins; i++) {
        new_tally->_batch_mu[i] *= amt;
        new_tally->_batch_variance[i] *= amt * amt;
        new_tally->_batch_std_dev[i] = sqrt(new_tally->_batch_variance[i]);
        new_tally->_batch_rel_err[i] = new_tally->_batch_std_dev[i] / 
                                        new_tally->_batch_mu[i];
    }

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally division with a constant operator.
 * @details This overloaded division operator allows the user to
 *          divide a tally by a constant. The creates a new DERIVED type
 *          tally, loads it with the dividend of the tally bin averages and the
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally / int(3)
 * @endcode
 *
 * @param amt the right operand and the constant value to divide the tally by
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator/(const int amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Divide the batch mu for this tally and updated the perturbed batch
     * statistics */
    for (int i=0; i < _num_bins; i++) {
        new_tally->_batch_mu[i] /= amt;
        new_tally->_batch_variance[i] *= (1.0 / amt) * (1.0 / amt);
        new_tally->_batch_std_dev[i] = sqrt(new_tally->_batch_variance[i]);
        new_tally->_batch_rel_err[i] = new_tally->_batch_std_dev[i] / 
                                        new_tally->_batch_mu[i];
    }

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally addition with a constant operator.
 * @details This overloaded addition operator allows the user to
 *          add a constant to a tally. The creates a new DERIVED type
 *          tally, loads it with the sum of the tally bin averages and 
 *          the constant, and updates its batch statistics appropriately. This
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally + float(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to add to the tally
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::operator+(const float amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    for (int i=0; i < _num_bins; i++) 
        new_tally->_batch_mu[i] += amt;

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally subtraction with a constant operator.
 * @details This overloaded subtraction operator allows the user to
 *          subtract a constant from a tally. The creates a new DERIVED type
 *          tally, loads it with the difference of the tally bin averages and 
 *          the constant, and updates its batch statistics appropriately. This
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally - float(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to subtract from the 
 *        tally
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::operator-(const float amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Subtracting an amount to this tally only subtracts the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    for (int i=0; i < _num_bins; i++) 
        new_tally->_batch_mu[i] -= amt;

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally multiplication with a constant operator.
 * @details This overloaded subtraction operator allows the user to
 *          multiply a constant with a tally. The creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally * float(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to multiply the tally by
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator*(const float amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Multiplying this tally by some amount multiplies the batch mu
     * by a constant and updates the perturbed batch statistics */
    for (int i=0; i < _num_bins; i++) {
        new_tally->_batch_mu[i] *= amt;
        new_tally->_batch_variance[i] *= amt * amt;
        new_tally->_batch_std_dev[i] = sqrt(new_tally->_batch_variance[i]);
        new_tally->_batch_rel_err[i] = new_tally->_batch_std_dev[i] / 
                                        new_tally->_batch_mu[i];
    }

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally division with a constant operator.
 * @details This overloaded division operator allows the user to
 *          divide a tally by a constant. The creates a new DERIVED type
 *          tally, loads it with the dividend of the tally bin averages and the
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally / float(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to divide the tally by
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator/(const float amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Divide the batch mu for this tally and updated the perturbed batch
     * statistics */
    for (int i=0; i < _num_bins; i++) {
        new_tally->_batch_mu[i] /= amt;
        new_tally->_batch_variance[i] *= (1.0 / amt) * (1.0 / amt);
        new_tally->_batch_std_dev[i] = sqrt(new_tally->_batch_variance[i]);
        new_tally->_batch_rel_err[i] = new_tally->_batch_std_dev[i] / 
                                        new_tally->_batch_mu[i];
    }

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally addition with a constant operator.
 * @details This overloaded addition operator allows the user to
 *          add a constant to a tally. The creates a new DERIVED type
 *          tally, loads it with the sum of the tally bin averages and 
 *          the constant, and updates its batch statistics appropriately. This
 *          method is intended to allow for simple tally arithmetic in Python.
 * 
 * @code 
 *          new_tally = tally + double(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to add to the tally
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::operator+(const double amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    for (int i=0; i < _num_bins; i++) 
        new_tally->_batch_mu[i] += amt;

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally subtraction with a constant operator.
 * @details This overloaded subtraction operator allows the user to
 *          subtract a constant from a tally. The creates a new DERIVED type
 *          tally, loads it with the difference of the tally bin averages and 
 *          the constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally - double(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to subtract from the 
 *        tally
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::operator-(const double amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Subtracting an amount to this tally only subtracts the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    for (int i=0; i < _num_bins; i++) 
        new_tally->_batch_mu[i] -= amt;

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally multiplication with a constant operator.
 * @details This overloaded subtraction operator allows the user to
 *          multiply a constant with a tally. The creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          constant, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python.
 *
 * @code 
 *          new_tally = tally * double(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to multiply the tally by
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator*(const double amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Multiplying this tally by some amount multiplies the batch mu
     * by a constant and updates the perturbed batch statistics */
    for (int i=0; i < _num_bins; i++) {
        new_tally->_batch_mu[i] *= amt;
        new_tally->_batch_variance[i] *= amt * amt;
        new_tally->_batch_std_dev[i] = sqrt(new_tally->_batch_variance[i]);
        new_tally->_batch_rel_err[i] = new_tally->_batch_std_dev[i] / 
                                        new_tally->_batch_mu[i];
    }

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally division with a constant operator.
 * @details This overloaded division operator allows the user to
 *          divide a tally by a constant. The creates a new DERIVED type
 *          tally, loads it with the dividend of the tally bin averages and the
 *          constant, and updates its batch statistics appropriately.
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *
 * @code 
 *          new_tally = tally / double(3.5)
 * @endcode
 *
 * @param amt the right operand and the constant value to divide the tally by
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::operator/(const double amt) {

    Tally* new_tally = this->clone();
    new_tally->_tally_domain = UNDEFINED;
    new_tally->_tally_type = DERIVED;   
    
    /* Divide the batch mu for this tally and updated the perturbed batch
     * statistics */
    for (int i=0; i < _num_bins; i++) {
        new_tally->_batch_mu[i] /= amt;
        new_tally->_batch_variance[i] *= (1.0 / amt) * (1.0 / amt);
        new_tally->_batch_std_dev[i] = sqrt(new_tally->_batch_variance[i]);
        new_tally->_batch_rel_err[i] = new_tally->_batch_std_dev[i] / 
                                        new_tally->_batch_mu[i];
    }

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally addition with an array of integers.
 * @details This overloaded division operator allows the user to
 *          ad a tally by an array of integers, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the sum of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          integer array and the length of the array - but in Python it only
 *          requires the integer array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2])
 *          new_tally = tally.addIntegers(array)
 * @endcode
 *
 * @param amt the right operand and the the array to add to the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::addIntegers(const int* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to add an integer array of length %d"
                        " to tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    double* batch_mu = new double[length];
    for (int i=0; i < _num_bins; i++)
        batch_mu[i] = _batch_mu[i] + amt[i];

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally addition with an array of floats.
 * @details This overloaded division operator allows the user to
 *          add a tally by an array of floats, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the sum of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          float array and the length of the array - but in Python it only
 *          requires the integer array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], numpy.dtype=float32)
 *          new_tally = tally.addFloats(array)
 * @endcode
 *
 * @param amt the right operand and the array to add to the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::addFloats(const float* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to add a float array of length %d"
                        " to tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    double* batch_mu = new double[length];
    for (int i=0; i < _num_bins; i++)
        batch_mu[i] = _batch_mu[i] + amt[i];

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);

}


/**
 * @brief Tally addition with an array of doubles.
 * @details This overloaded division operator allows the user to
 *          add a tally to an array of doubles, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the sum of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          double array and the length of the array - but in Python it only
 *          requires the double array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], numpy.dtype=float64)
 *          new_tally = tally.addDoubles(array)
 * @endcode
 *
 * @param amt the right operand and the array to add to the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally sum
 */
DerivedTally* Tally::addDoubles(const double* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to add a float array of length %d"
                        " to tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    double* batch_mu = new double[length];
    for (int i=0; i < _num_bins; i++)
        batch_mu[i] = _batch_mu[i] + amt[i];

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);

}


/**
 * @brief Tally subtraction with an array of integers.
 * @details This overloaded subtraction operator allows the user to
 *          subtract an array of integers from a tally, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the difference of the tally bin averages and
 *           the array, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          integer array and the length of the array - but in Python it only
 *          requires the integer array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.integer)
 *          new_tally = tally.subtractIntegers(array)
 * @endcode
 *
 * @param amt the right operand and the array to subtract from the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::subtractIntegers(const int* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to subtract an integer array of length %d"
                        " to tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Subtracting an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    double* batch_mu = new double[length];
    for (int i=0; i < _num_bins; i++)
        batch_mu[i] = _batch_mu[i] - amt[i];

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally subtraction with an array of floats.
 * @details This overloaded subtraction operator allows the user to
 *          subtract an array of floats from a tally, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the difference of the tally bin averages and 
 *          the array, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          float array and the length of the array - but in Python it only
 *          requires the float array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.float32)
 *          new_tally = tally.subtractFloats(array)
 * @endcode
 *
 * @param amt the right operand and the array to subtract from the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::subtractFloats(const float* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to subtract a float array of length %d"
                        " to tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    double* batch_mu = new double[length];
    for (int i=0; i < _num_bins; i++)
        batch_mu[i] = _batch_mu[i] - amt[i];

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally subtraction with an array of doubles.
 * @details This overloaded subtraction operator allows the user to
 *          subtract an array of doubles from a tally, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the difference of the tally bin averages and 
 *          the array, and updates its batch statistics appropriately. This 
 *          method is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          double array and the length of the array - but in Python it only
 *          requires the double array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.float64)
 *          new_tally = tally.subtractIntegers(array)
 * @endcode
 *
 * @param amt the right operand and the array to subtract from the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally difference
 */
DerivedTally* Tally::subtractDoubles(const double* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to subtract a double array of length %d"
                        " to tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and does not otherwise perturb the batch statistics */
    double* batch_mu = new double[length];
    for (int i=0; i < _num_bins; i++)
        batch_mu[i] = _batch_mu[i] - amt[i];

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally multiplication with an array of integers.
 * @details This overloaded multiplication operator allows the user to
 *          multiply an array of integers with a tally, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          integer array and the length of the array - but in Python it only
 *          requires the integer array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.int32)
 *          new_tally = tally.multiplyIntegers(array)
 * @endcode
 *
 * @param amt the right operand and the array to multiply with the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::multiplyIntegers(const int* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to multiply an integer array of length %d"
                        " with tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and adjust batch statistics appropriately */
    double* batch_mu = new double[length];
    double* batch_variance = new double[length];
    double* batch_std_dev = new double[length];
    double* batch_rel_err = new double[length];

    for (int i=0; i < _num_bins; i++) {
        batch_mu[i] = _batch_mu[i] * amt[i];
        batch_variance[i] = amt[i] * amt[i] * _batch_variance[i];
        batch_std_dev[i] = sqrt(batch_variance[i]);
        batch_rel_err[i] = batch_std_dev[i] / batch_mu[i];
    }

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setBatchVariance(batch_variance);
    static_cast<DerivedTally*>(new_tally)->setBatchStdDev(batch_std_dev);
    static_cast<DerivedTally*>(new_tally)->setBatchRelErr(batch_rel_err);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally multiplication with an array of floats.
 * @details This overloaded multiplication operator allows the user to
 *          multiply an array of floats with a tally, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          float array and the length of the array - but in Python it only
 *          requires the float array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.float32)
 *          new_tally = tally.multiplyFloats(array)
 * @endcode
 *
 * @param amt the right operand and the array to multiply with the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::multiplyFloats(const float* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to multiply a float array of length %d"
                        " with tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and adjust batch statistics appropriately */
    double* batch_mu = new double[length];
    double* batch_variance = new double[length];
    double* batch_std_dev = new double[length];
    double* batch_rel_err = new double[length];

    for (int i=0; i < _num_bins; i++) {
        batch_mu[i] = _batch_mu[i] * amt[i];
        batch_variance[i] = amt[i] * amt[i] * _batch_variance[i];
        batch_std_dev[i] = sqrt(batch_variance[i]);
        batch_rel_err[i] = batch_std_dev[i] / batch_mu[i];
    }

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setBatchVariance(batch_variance);
    static_cast<DerivedTally*>(new_tally)->setBatchStdDev(batch_std_dev);
    static_cast<DerivedTally*>(new_tally)->setBatchRelErr(batch_rel_err);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally multiplication with an array of doubles.
 * @details This overloaded multiplication operator allows the user to
 *          multiply an array of doubles with a tally, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          double array and the length of the array - but in Python it only
 *          requires the double array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.float64)
 *          new_tally = tally.multiplyDoubles(array)
 * @endcode
 *
 * @param amt the right operand and the array to multiply with the tally
 * @param length the length of the array
 * @return a DERIVED type tally with the tally product
 */
DerivedTally* Tally::multiplyDoubles(const double* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to multiply a double array of length %d"
                        " with tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);

    /* Adding an amount to this tally only increments the batch mu
     * by an offset and adjust batch statistics appropriately */
    double* batch_mu = new double[length];
    double* batch_variance = new double[length];
    double* batch_std_dev = new double[length];
    double* batch_rel_err = new double[length];

    for (int i=0; i < _num_bins; i++) {
        batch_mu[i] = _batch_mu[i] * amt[i];
        batch_variance[i] = amt[i] * amt[i] * _batch_variance[i];
        batch_std_dev[i] = sqrt(batch_variance[i]);
        batch_rel_err[i] = batch_std_dev[i] / batch_mu[i];
    }


    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setBatchVariance(batch_variance);
    static_cast<DerivedTally*>(new_tally)->setBatchStdDev(batch_std_dev);
    static_cast<DerivedTally*>(new_tally)->setBatchRelErr(batch_rel_err);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally division with an array of integers.
 * @details This overloaded multiplication operator allows the user to
 *          divide a tally by an array of integers, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          integer array and the length of the array - but in Python it only
 *          requires the integer array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.int32)
 *          new_tally = tally.divideIntegers(array)
 * @endcode
 *
 * @param amt the right operand and the array to divide the tally by
 * @param length the length of the array
 * @return a DERIVED type tally with the tally dividend
 */
DerivedTally* Tally::divideIntegers(const int* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to divide an integer array of length %d"
                        " with tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and adjust batch statistics appropriately */
    double* batch_mu = new double[length];
    double* batch_variance = new double[length];
    double* batch_std_dev = new double[length];
    double* batch_rel_err = new double[length];

    for (int i=0; i < _num_bins; i++) {
        batch_mu[i] = (_batch_mu[i] / amt[i]);
        batch_variance[i] = (1.0 / amt[i]) * (1.0 / amt[i]) 
                                                        * _batch_variance[i];
        batch_std_dev[i] = sqrt(batch_variance[i]);
        batch_rel_err[i] = batch_std_dev[i] / batch_mu[i];
    }

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setBatchVariance(batch_variance);
    static_cast<DerivedTally*>(new_tally)->setBatchStdDev(batch_std_dev);
    static_cast<DerivedTally*>(new_tally)->setBatchRelErr(batch_rel_err);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally division with an array of floats.
 * @details This overloaded multiplication operator allows the user to
 *          divide a tally by an array of floats, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          float array and the length of the array - but in Python it only
 *          requires the float array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.float32)
 *          new_tally = tally.divideFloats(array)
 * @endcode
 *
 * @param amt the right operand and the array to divide the tally by
 * @param length the length of the array
 * @return a DERIVED type tally with the tally dividend
 */
DerivedTally* Tally::divideFloats(const float* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to divide a float array of length %d"
                        " with tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and adjust batch statistics appropriately */
    double* batch_mu = new double[length];
    double* batch_variance = new double[length];
    double* batch_std_dev = new double[length];
    double* batch_rel_err = new double[length];

    for (int i=0; i < _num_bins; i++) {
        batch_mu[i] = (_batch_mu[i] / amt[i]);
        batch_variance[i] = (1.0 / amt[i]) * (1.0 / amt[i]) 
                                                        * _batch_variance[i];
        batch_std_dev[i] = sqrt(batch_variance[i]);
        batch_rel_err[i] = batch_std_dev[i] / batch_mu[i];
    }

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setBatchVariance(batch_variance);
    static_cast<DerivedTally*>(new_tally)->setBatchStdDev(batch_std_dev);
    static_cast<DerivedTally*>(new_tally)->setBatchRelErr(batch_rel_err);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/**
 * @brief Tally division with an array of doubles.
 * @details This overloaded multiplication operator allows the user to
 *          divide a tally by an array of doubles, if the number of values
 *          equals the number of tally bins. This creates a new DERIVED type
 *          tally, loads it with the product of the tally bin averages and the
 *          array, and updates its batch statistics appropriately. This method
 *          is intended to allow for simple tally arithmetic in Python. 
 *           This method prototype appears to require two operands - the 
 *          double array and the length of the array - but in Python it only
 *          requires the double array as follows:
 * 
 * @code 
 *          array = numpy.array([1.2, 3.5, 4.7, 8.2], dtype=numpy.float64)
 *          new_tally = tally.divideDoubles(array)
 * @endcode
 *
 * @param amt the right operand and the array to divide the tally by
 * @param length the length of the array
 * @return a DERIVED type tally with the tally dividend
 */
DerivedTally* Tally::divideDoubles(const double* amt, const int length) {

    if (length != _num_bins)
        log_printf(ERROR, "Unable to divide a double array of length %d"
                        " with tally %s with %d bins", length,
                        _tally_name, _num_bins);

    Tally* new_tally = this->clone();
    new_tally->setTallyDomainType(UNDEFINED);
    new_tally->setTallyType(DERIVED);
    new_tally->setNumBatches(1);
    
    /* Adding an amount to this tally only increments the batch mu
     * by an offset and adjust batch statistics appropriately */
    double* batch_mu = new double[length];
    double* batch_variance = new double[length];
    double* batch_std_dev = new double[length];
    double* batch_rel_err = new double[length];

    for (int i=0; i < _num_bins; i++) {
        batch_mu[i] = (_batch_mu[i] / amt[i]);
        batch_variance[i] = (1.0 / amt[i]) * (1.0 / amt[i]) 
                                                        * _batch_variance[i];
        batch_std_dev[i] = sqrt(batch_variance[i]);
        batch_rel_err[i] = batch_std_dev[i] / batch_mu[i];
    }

    new_tally->setNumBatches(1);
    static_cast<DerivedTally*>(new_tally)->setBatchMu(batch_mu);
    static_cast<DerivedTally*>(new_tally)->setBatchVariance(batch_variance);
    static_cast<DerivedTally*>(new_tally)->setBatchStdDev(batch_std_dev);
    static_cast<DerivedTally*>(new_tally)->setBatchRelErr(batch_rel_err);
    static_cast<DerivedTally*>(new_tally)->setComputedBatchStatistics(true);

    return static_cast<DerivedTally*>(new_tally);
}


/*****************************************************************************/
/********************** Collision Rate Tallying Methods **********************/
/*****************************************************************************/

/**
 * @brief Tally the isotope collision rate by incrementing the tally by 1.0
 *        at the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeCollisionRateTally::tally(neutron* neutron) {
    Tally::tally(neutron, 1.0);
    return;
}


/**
 * @brief Tally the material collision rate by incrementing the tally by 1.0
 *         at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialCollisionRateTally::tally(neutron* neutron) {
    Tally::tally(neutron, 1.0);
    return;
}


/**
 * @brief Tally the region collision rate by incrementing the tally by 1.0
 *         at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionCollisionRateTally::tally(neutron* neutron) {
    Tally::tally(neutron, 1.0);
    return;
}


/**
 * @brief Tally the geometry collision rate by incrementing the tally by 1.0
 *         at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryCollisionRateTally::tally(neutron* neutron) {
    Tally::tally(neutron, 1.0);
    return;
}

/******************************************************************************/
/************************ Elastic Rate Tallying Methods ***********************/
/******************************************************************************/
//FIXME

/**
 * @brief Tally the isotope elastic scattering rate by incrementing the 
 *        tally by \f$ \frac{\sigma_s}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeElasticRateTally::tally(neutron* neutron) {
    double weight = _isotope->getElasticXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material elastic scattering rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_s}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialElasticRateTally::tally(neutron* neutron) {
    double weight = _material->getElasticMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material elastic scattering rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_s}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionElasticRateTally::tally(neutron* neutron) {
    double weight = _region->getElasticMacroXS(neutron->_old_energy) 
                   / neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the geometry elastic scattering rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_s}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryElasticRateTally::tally(neutron* neutron) {
    double weight = neutron->_region->getElasticMacroXS(neutron->_old_energy)
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}

/******************************************************************************/
/********************** Absorption Rate Tallying Methods **********************/
/******************************************************************************/

/**
 * @brief Tally the isotope absorption rate by incrementing the 
 *        tally by \f$ \frac{\sigma_a}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeAbsorptionRateTally::tally(neutron* neutron) {
    double weight = _isotope->getAbsorptionXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material absorption rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_a}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialAbsorptionRateTally::tally(neutron* neutron) {
    double weight = _material->getAbsorptionMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the region absorption rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_a}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionAbsorptionRateTally::tally(neutron* neutron) {
    double weight = _region->getAbsorptionMacroXS(neutron->_old_energy) 
   		    / neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the region absorption rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_a}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryAbsorptionRateTally::tally(neutron* neutron) {
    double weight = neutron->_region->getAbsorptionMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}

/******************************************************************************/
/************************ Capture Rate Tallying Methods ***********************/
/******************************************************************************/

/**
 * @brief Tally the isotope capture rate by incrementing the 
 *        tally by \f$ \frac{\sigma_c}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeCaptureRateTally::tally(neutron* neutron) {
    double weight = _isotope->getCaptureXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material capture rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_c}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialCaptureRateTally::tally(neutron* neutron) {
    double weight = _material->getCaptureMacroXS(neutron->_old_energy) 
 			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the region capture rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_c}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionCaptureRateTally::tally(neutron* neutron) {
    double weight = _region->getCaptureMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the geometry capture rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_c}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryCaptureRateTally::tally(neutron* neutron) {
    double weight = neutron->_region->getCaptureMacroXS(neutron->_old_energy) 
  			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/******************************************************************************/
/*********************** Fission Rate Tallying Methods ************************/
/******************************************************************************/

/**
 * @brief Tally the isotope fission rate by incrementing the 
 *        tally by \f$ \frac{\sigma_f}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeFissionRateTally::tally(neutron* neutron) {
    double weight = _isotope->getFissionXS(neutron->_old_energy) 
 		     / neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material fission rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_f}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialFissionRateTally::tally(neutron* neutron) {
    double weight = _material->getFissionMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}

/**
 * @brief Tally the region fission rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_f}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionFissionRateTally::tally(neutron* neutron) {
    double weight = _region->getFissionMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the geometry fission rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_f}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryFissionRateTally::tally(neutron* neutron) {
    double weight = neutron->_region->getFissionMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/******************************************************************************/
/********************** Transport Rate Tallying Methods ***********************/
/******************************************************************************/

/**
 * @brief Tally the isotope transport rate by incrementing the 
 *        tally by \f$ \frac{\sigma_tr}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeTransportRateTally::tally(neutron* neutron) {
    double weight = _isotope->getTransportXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material transport rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_tr}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialTransportRateTally::tally(neutron* neutron) {
    double weight = _material->getTransportMacroXS(neutron->_old_energy) 
		     /neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the region transport rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_tr}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionTransportRateTally::tally(neutron* neutron) {
    double weight = _region->getTransportMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the geometry transport rate by incrementing the 
 *        tally by \f$ \frac{\Sigma_tr}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryTransportRateTally::tally(neutron* neutron) {
    double weight = neutron->_region->getTransportMacroXS(neutron->_old_energy) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}

/******************************************************************************/
/********************** Diffusion Rate Tallying Methods ***********************/
/******************************************************************************/

/**
 * @brief Tally the isotope diffusion rate by incrementing the 
 *        tally by \f$ \frac{\frac{1}{3\sigma_tr}}{\Sigma_t} \f$ at 
 *        the neutron's energy.
 * @param neutron the neutron of interest
 */
void IsotopeDiffusionRateTally::tally(neutron* neutron) {
    double weight = 1.0 / (3.0 * _isotope->getTransportXS(neutron->_old_energy))  		         / neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the material diffusion rate by incrementing the 
 *        tally by \f$ \frac{\frac{1}{3\Sigma_tr}}{\Sigma_t} \f$ at 
 *        the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialDiffusionRateTally::tally(neutron* neutron) {
    double weight = 1.0 / (3.0 * 
			   _material->getTransportMacroXS(neutron->_old_energy)) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the region diffusion rate by incrementing the 
 *        tally by \f$ \frac{\frac{1}{3\Sigma_tr}}{\Sigma_t} \f$ at 
 *        the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionDiffusionRateTally::tally(neutron* neutron) {
    double weight = 1.0 / (3.0 * 
			   _region->getTransportMacroXS(neutron->_old_energy)) 
			/ neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the geometry diffusion rate by incrementing the 
 *        tally by \f$ \frac{\frac{1}{3\Sigma_tr}}{\Sigma_t} \f$ at 
 *        the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryDiffusionRateTally::tally(neutron* neutron) {
    double weight = 1.0 / (3.0 * 
		   neutron->_region->getTransportMacroXS(neutron->_old_energy))
                  / neutron->_total_xs;
    Tally::tally(neutron, weight);
    return;
}


/******************************************************************************/
/*********************** Leakage Rate Tallying Methods ************************/
/******************************************************************************/
//FIXME

/**
 * @brief Tally the material leakage rate by incrementing the 
 *        tally by \f$ DB^2 \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialLeakageRateTally::tally(neutron* neutron) {
    double weight = _material->getBucklingSquared() / 
   		(3.0 * _material->getTransportMacroXS(neutron->_old_energy) * 
		neutron->_total_xs);
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the region leakage rate by incrementing the 
 *        tally by \f$ DB^2 \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionLeakageRateTally::tally(neutron* neutron) {
    double weight = _region->getBucklingSquared() / 
		(3.0 * _region->getTransportMacroXS(neutron->_old_energy) * 
		neutron->_total_xs);
    Tally::tally(neutron, weight);
    return;
}


/**
 * @brief Tally the geometry leakage rate by incrementing the 
 *        tally by \f$ DB^2 \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryLeakageRateTally::tally(neutron* neutron) {
    double weight = neutron->_region->getBucklingSquared() / 
	(3.0 * neutron->_region->getTransportMacroXS(neutron->_old_energy) * 
	neutron->_total_xs);
    Tally::tally(neutron, weight);
    return;
}


/******************************************************************************/
/************************** Flux Tallying Methods *****************************/
/******************************************************************************/

/**
 * @brief Tally the material flux by incrementing the tally 
 *        by \f$ \frac{1}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void MaterialFluxTally::tally(neutron* neutron) {
    Tally::tally(neutron, double(1.0 / neutron->_total_xs));
    return;
}


/**
 * @brief Tally the region flux by incrementing the tally 
 *        by \f$ \frac{1}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void RegionFluxTally::tally(neutron* neutron) {
     Tally::tally(neutron, double(1.0 / neutron->_total_xs));
    return;
}


/**
 * @brief Tally the geomety flux by incrementing the 
 *        tally by \f$ \frac{1}{\Sigma_t} \f$ at the neutron's energy.
 * @param neutron the neutron of interest
 */
void GeometryFluxTally::tally(neutron* neutron) {
    Tally::tally(neutron, double(1.0 / neutron->_total_xs));
    return;
}


/******************************************************************************/
/*************************** Mean Lifetime Tallies ****************************/
/******************************************************************************/

/**
 * @brief Tally the time between collisions by incrementing the tally with the
 *        travel time spent by this neutron to reach this collision:
 *
 *        \f$ d = \frac{1}{\Sigma_t} \f$
 *        \f$ v = c\sqrt{2Em_n} \f$
 *        \f$ time = \frac{d}{v} \f$
 * @param neutron the neutron of interest
 */
void MaterialInterCollisionTimeTally::tally(neutron* neutron) {
    float distance = (1.0 / neutron->_total_xs) * 1E-2;
    float velocity = LIGHT_SPEED * sqrt(2.0 * neutron->_old_energy 
					/ NEUTRON_MASS);
    Tally::tally(neutron, double(distance / velocity));
    return;
}



/**
 * @brief Tally the time between collisions by incrementing the tally with the
 *        travel time spent by this neutron to reach this collision:
 *
 *        \f$ d = \frac{1}{\Sigma_t} \f$
 *        \f$ v = c\sqrt{2Em_n} \f$
 *        \f$ time = \frac{d}{v} \f$
 * @param neutron the neutron of interest
 */
void RegionInterCollisionTimeTally::tally(neutron* neutron) {
    float distance = (1.0 / neutron->_total_xs)  * 1E-2;
    float velocity = LIGHT_SPEED * sqrt(2.0 * neutron->_old_energy 
					/ NEUTRON_MASS);
    Tally::tally(neutron, double(distance / velocity));
    return;
}


/**
 * @brief Tally the time between collisions by incrementing the tally with the
 *        travel time spent by this neutron to reach this collision:
 *
 *        \f$ d = \frac{1}{\Sigma_t} \f$
 *        \f$ v = c\sqrt{2Em_n} \f$
 *        \f$ time = \frac{d}{v} \f$
 * @param neutron the neutron of interest
 */
void GeometryInterCollisionTimeTally::tally(neutron* neutron) {
    float distance = (1.0 / neutron->_total_xs) * 1E-2;
    float velocity = LIGHT_SPEED * sqrt(2.0 * neutron->_old_energy 
					/ NEUTRON_MASS);
    Tally::tally(neutron, double(distance / velocity));
    return;
}


/******************************************************************************/
/****************************** Derived Tallies *******************************/
/******************************************************************************/

/**
 * @brief A DERIVED tally cannot tally anything, and will throw an exception.
 * @param neutron the neutron of interest
 */
void DerivedTally::tally(neutron* neutron) {
    log_printf(ERROR, "Unable to tally a neutron in DERIVED type tally %s", 
	       _tally_name);
}


/**
 * @brief Sets the tally name.
 * @param tally_name the name of the tally.
 */
void DerivedTally::setTallyName(char* tally_name) {
    _tally_name = tally_name;
}


/**
 * @brief Assigns an array for the tally values.
 * @param tallies an array of tally values
 */
void DerivedTally::setTallies(double** tallies) {
    _tallies = tallies;
}


/**
 * @brief Assigns an array for the tally batch averages.
 * @param batch_mu an array of tally averages
 */
void DerivedTally::setBatchMu(double* batch_mu) {
    _batch_mu = batch_mu;
}


/**
 * @brief Assigns an array for the tally batch variances.
 * @param batch_variance an array of tally variances
 */
void DerivedTally::setBatchVariance(double* batch_variance) {
    _batch_variance = batch_variance;
}



/**
 * @brief Assigns an array for the tally batch standard deviations.
 * @param batch_std_dev an array of tally standard deviations
 */
void DerivedTally::setBatchStdDev(double* batch_std_dev) {
    _batch_std_dev = batch_std_dev;
}



/**
 * @brief Assigns an array for the tally batch relative errors.
 * @param batch_rel_err an array of tally relative errors
 */
void DerivedTally::setBatchRelErr(double* batch_rel_err) {
    _batch_rel_err = batch_rel_err;
}



/**
 * @brief Assigns whether or not the tally has computed batch statistics.
 * @param computed true if batch statistics are computed; otherwise false
 */
void DerivedTally::setComputedBatchStatistics(bool computed) {
    _computed_statistics = computed;
}
