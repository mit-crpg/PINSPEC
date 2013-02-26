/*
 * tally.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TALLY_H_
#define TALLY_H_

#include <limits>
#include <math.h>
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
	REGION
} tallyDomainType;


/* Type of tallies */
typedef enum tallyTypes {
	COLLISION,
	FLUX,
	ELASTIC,
	ABSORPTION,
	CAPTURE,
	FISSION,
	TRANSPORT,
	DIFFUSION,
	LEAKAGE
} tallyType;


/**
 * This class represents a set of tallies. A set of values
 * define the edges between bins for each tally. This class 
 * holds the edges, the centers between bins. It also allows 
 * for tallies to be made within each bin.
 */
class Tally{

private:
	char* _name;
	int _num_bins;
	float* _edges;
	double* _centers;
	double* _tallies;
	int* _num_tallies;
	float _bin_delta;
	binSpacingType _bin_spacing;
	tallyDomainType _tally_domain;
	tallyType _tally_type;
	char* _isotopes;

public:
	Tally();
	virtual ~Tally();
	char* getTallyName();
	int getNumBins();
	float* getBinEdges();
	double* getBinCenters();
	float getBinDelta();
	float getBinDelta(float sample);
	binSpacingType getBinSpacingType();
	tallyDomainType getTallyDomainType();
	tallyType getTallyType();
	double* getTallies();
	double getTally(int bin_index);
	int* getNumTallies();
	int getNumTallies(int bin_index);
	double getMaxTally();
	double getMinTally();
	int getBinIndex(float sample);
	char* getIsotopes();

	void setTallyName(char* name);
	void setTallyDomainType(tallyDomainType type);
	void setTallyType(tallyType type);
	void setBinEdges(float* edges, int num_edges);
	void setIsotopes(char* isotopes);

	void generateBinEdges(float start, float end, int num_bins, 
											binSpacingType type);
	void generateBinCenters();
	void tally(float* samples, int num_samples);
	void tally(float sample);
	void weightedTally(float* samples, float* sample_weights, int num_samples);
	void weightedTally(float sample, float weight);
	void normalizeTallies();
	void normalizeTallies(float scale_factor);
};

#endif /* BINNER_H_ */
