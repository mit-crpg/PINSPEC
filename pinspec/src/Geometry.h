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

#ifdef __cplusplus
#include <sys/stat.h>
#include <omp.h>
#include <vector>
#include "Region.h"
#include "Fissioner.h"
#include "TallyBank.h"
#include "Timer.h"


typedef enum spatialTypes {
	INFINITE_HOMOGENEOUS,
	HOMOGENEOUS_EQUIVALENCE,
	HETEROGENEOUS,
} spatialType;


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

	float _dancoff;
	float _sigma_e;
	float _beta;
	float _alpha1;
	float _alpha2;

	float _buckling_squared;

    /* Ratio of sigam_f * vol_f / sigma_m * vol_m */
    /* Compute ratios ahead of time as an optimization */
    int _num_ratios;
    float* _pmf_ratios;
	binSpacingTypes _scale_type;
	float _start_energy;
	float _end_energy;
	float _delta_energy;

    void initializeProbModFuelRatios();
	int getEnergyGridIndex(float energy) const;
    float computeFuelFuelCollisionProb(neutron* neutron);
    float computeModeratorFuelCollisionProb(neutron* neutron);

public:
	Geometry();
	virtual ~Geometry();
	
	/* getters */
	int getNumNeutronsPerBatch();
	int getTotalNumNeutrons();
	int getNumBatches();
	int getNumThreads();
	spatialType getSpatialType();
	float getBucklingSquared();
			
	/* setters */
	void setNeutronsPerBatch(int num_neutrons_per_batch);
	void setNumBatches(int num_batches);
	void setNumThreads(int num_threads);
	void setSpatialType(spatialType spatial_type);
	void setDancoffFactor(float dancoff);
	void addRegion(Region* region);
	void setBucklingSquared(float buckling_squared);

	/* Monte Carlo kernel */
	void runMonteCarloSimulation();
};


/**
 * This method returns the index for a certain energy (eV) into
 * the uniform lethargy grid of the Geometry's Pmf ratios
 * @param energy the energy (eV) of interest
 * @return the index into the uniform lethargy grid
 */
inline int Geometry::getEnergyGridIndex(float energy) const {

	int index;

	energy = log10(energy);

	if (energy > _end_energy)
		index = _num_ratios - 1;
	else if (energy < _start_energy)
		index = 0;
	else
		index = int(floor((energy - _start_energy) / _delta_energy));

	return index;
}

#endif

#endif /* GEOMETRY_H_ */
