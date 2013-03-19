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
#include <omp.h>
#include <vector>
#include "Region.h"
#include "Fissioner.h"
#include "Neutron.h"
#include "Tally.h"
#include "Timer.h"       


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

    float _dancoff;
	float _sigma_e;
	float _beta;
	float _alpha1;
	float _alpha2;

    /* Ratio of sigam_f * vol_f / sigma_m * vol_m */
    /* Compute ratios ahead of time as an optimization */
    int _num_ratios;
    float* _pmf_ratios;
	binSpacingTypes _scale_type;
	float _start_energy;
	float _end_energy;
	float _delta_energy;

    std::vector<Tally*> _tallies;

	float getTotalMacroXS(float energy, Region* region);
	float getTotalMacroXS(int energy_index, Region* region);
	float getTotalMicroXS(float energy, Region* region);
	float getTotalMicroXS(int energy_index, Region* region);
	
	float getElasticMacroXS(float energy, Region* region);
	float getElasticMacroXS(int energy_index, Region* region);
	float getElasticMicroXS(float energy, Region* region);
	float getElasticMicroXS(int energy_index, Region* region);
	
	float getAbsorptionMacroXS(float energy, Region* region);
	float getAbsorptionMacroXS(int energy_index, Region* region);
	float getAbsorptionMicroXS(float energy, Region* region);
	float getAbsorptionMicroXS(int energy_index, Region* region);
	
	float getCaptureMacroXS(float energy, Region* region);
	float getCaptureMacroXS(int energy_index, Region* region);
	float getCaptureMicroXS(float energy, Region* region);
	float getCaptureMicroXS(int energy_index, Region* region);
	
	float getFissionMacroXS(float energy, Region* region);
	float getFissionMacroXS(int energy_index, Region* region);
	float getFissionMicroXS(float energy, Region* region);
	float getFissionMicroXS(int energy_index, Region* region);
	
	float getTransportMicroXS(float energy, Region* region);
	float getTransportMicroXS(int energy_index, Region* region);
	float getTransportMacroXS(float energy, Region* region);
	float getTransportMacroXS(int energy_index, Region* region);

    void tally(float sample, int batch_num, Region* region, collisionType type);
    void initializeBatchTallies();
    void initializeProbModFuelRatios();
    bool isPrecisionTriggered();
    void incrementNumBatches(int num_batches);
	int getEnergyGridIndex(float energy) const;
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
    void setDancoffFactor(float dancoff);
	void addRegion(Region* region);
    void addTally(Tally* tally);

	/* Monte Carlo kernel */
	void runMonteCarloSimulation();
    void computeBatchStatistics();
    void computeScaledBatchStatistics();
    void outputBatchStatistics(char* directory,  char* suffix);
};


/**
 * This method returns the index for a certain energy (eV) into
 * the uniform energy grid of the Geometry's Pmf ratios
 * @param energy the energy (eV) of interest
 * @return the index into the uniform energy grid
 */
inline int Geometry::getEnergyGridIndex(float energy) const {

	int index;

	if (_spatial_type != HOMOGENEOUS_EQUIVALENCE)
	    log_printf(ERROR, "Unable to return an index for pmf ratios "
		       			"since the geometry is not HOMOGENEOUS_EQUIVALENCE"
						" spatial type");

	if (_scale_type == EQUAL) {
		if (energy > _end_energy)
			index = _num_ratios - 1;
		else if (energy < _start_energy)
			index = 0;
		else
			index = floor((energy - _start_energy) / _delta_energy);
	}

	else if (_scale_type == LOGARITHMIC)
		energy = log10(energy);

		if (energy > _end_energy)
			index = _num_ratios - 1;
		else if (energy < _start_energy)
			index = 0;
		else
			index = floor((energy - _start_energy) / _delta_energy);

	return index;
}

#endif

#endif /* GEOMETRY_H_ */
