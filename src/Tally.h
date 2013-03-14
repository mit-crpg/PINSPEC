/*
 * Tally.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TALLY_H_
#define TALLY_H_

#ifdef __cplusplus
#include <limits>
#include <math.h>
#endif
#include "log.h"
#include "arraycreator.h"


/* Tally spacing types */
typedef enum binSpacingTypes {
	EQUAL,
	LOGARITHMIC,
	OTHER
} binSpacingType;


/* The domain in which this tally resides */
typedef enum tallyDomainTypes {
	MATERIAL,
	ISOTOPE,
	REGION,
    GEOMETRY
} tallyDomainType;


/* Type of tallies */
typedef enum tallyTypes {
	FLUX,
	COLLISION_RATE,
	ELASTIC_RATE,
	ABSORPTION_RATE,
	CAPTURE_RATE,
	FISSION_RATE,
	TRANSPORT_RATE,
	DIFFUSION_RATE,
	LEAKAGE_RATE
} tallyType;


/**
 * This class represents a set of tallies. A set of values
 * define the edges between bins for each tally. This class 
 * holds the edges, the centers between bins. It also allows 
 * for tallies to be made within each bin.
 */
#ifdef __cplusplus
class Tally{

private:
	char* _tally_name;
	int _num_bins;
	double* _edges;
	double* _centers;
	double** _tallies;
	int** _num_tallies;
	double _bin_delta;
	binSpacingType _bin_spacing;
	tallyDomainType _tally_domain;
	tallyType _tally_type;

	int _num_batches;
	double* _batch_mu;
	double* _batch_variance;
	double* _batch_std_dev;
	double* _batch_rel_err;
	bool _computed_statistics;

public:
	Tally(char* tally_name, tallyDomainType tally_domain, tallyType tally_type);
	virtual ~Tally();
	char* getTallyName();
	int getNumBins();
	double* getBinEdges();
	double* getBinCenters();
	double getBinDelta();
	double getBinDelta(double sample);
	binSpacingType getBinSpacingType();
	tallyDomainType getTallyDomainType();
	tallyType getTallyType();
	double** getTallies();
	double getTally(int bin_index, int batch_num);
	int** getNumTallies();
	int getNumTallies(int batch_num, int bin_index);
	double getMaxTally();
	double getMinTally();
	int getBinIndex(double sample);

    /* IMPORTANT: The following five class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be transformed
     * into numpy arrays */
    void retrieveTallyEdges(double* data, int num_bins);
    void retrieveTallyCenters(double* data, int num_bins);
    void retrieveTallyMu(double* data, int num_bins);
    void retrieveTallyVariance(double* data, int num_bins);
    void retrieveTallyStdDev(double* data, int num_bins);
    void retrieveTallyRelErr(double* data, int num_bins);

	int getNumBatches();
	double* getBatchMu();
	double* getBatchVariance();
	double* getBatchStdDev();
	double* getBatchRelativeError();

    void setBinSpacingType(binSpacingType type);
	void setBinEdges(double* edges, int num_edges);

	void generateBinEdges(double start, double end, int num_bins, 
											binSpacingType type);
	void setNumBatches(int num_batches);
	Tally* clone();

	void generateBinCenters();

	void tally(double* samples, int num_samples, int batch_num);
	void tally(double sample, int batch_num);
	void weightedTally(double* samples, double* sample_weights, 
			   int num_samples, int batch_num);
	void weightedTally(double sample, double weight, int batch_num);
	void normalizeTallies();
	void normalizeTallies(double scale_factor);
	void computeBatchStatistics();
	void computeScaledBatchStatistics(double scale_factor);
	void outputBatchStatistics(const char* filename);
};

#endif

#endif /* BINNER_H_ */
