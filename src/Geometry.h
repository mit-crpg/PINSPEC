/*
 * Geometry.h
 *
 *  Created on: Mar 6, 2013
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */	

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <sys/stat.h>
//#include <omp.h>
#include "/Users/samuelshaner/Downloads/gcc-4.6.2/gcc-4.6.2/build/prev-x86_64-apple-darwin11.2.0/i386/libgomp/omp.h"
#include "Region.h"
#include "Fissioner.h"
#include "Neutron.h"


typedef enum spatialTypes {
	INFINITE_HOMOGENEOUS,
	HOMOGENEOUS_EQUIVALENCE,
	HETEROGENEOUS,
} spatialType;


#ifdef __cplusplus
class Geometry {

private:

	int _num_neutrons_per_batch;
	int _num_batches;
	int _num_threads;

	spatialType _spatial_type;
	Region* _infinite_medium;
	Region* _fuel;
	Region* _moderator;
	neutron* _neutrons;
	Fissioner* _fissioner;

	float _sigma_e;
	float _beta;
	float _alpha1;
	float _alpha2;

    float computeFuelFuelCollisionProb(float energy);
    float computeModeratorFuelCollisionProb(float energy);

public:
	Geometry();
	virtual ~Geometry();
	
	/* getters */
	int getNumNeutronsPerBatch();
	int getTotalNumNeutrons();
	int getNumBatches();
	int getNumThreads();
	spatialType getSpatialType();
			
	/* setters */
	void setNeutronsPerBatch(int num_neutrons_per_batch);
	void setNumBatches(int num_batches);
	void setNumThreads(int num_threads);
	void setSpatialType(spatialType spatial_type);
	void setTwoRegionPinCellParams(float sigma_e, float beta,
								float alpha1, float alpha2);
	void addRegion(Region* region);

	/* Monte Carlo kernel */
	void runMonteCarloSimulation();
    void computeBatchStatistics();
    void computeScaledBatchStatistics();
    void outputBatchStatistics(char* directory,  char* suffix);
};

#endif

#endif /* GEOMETRY_H_ */
