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

//#include <vector>
//#include "log.h"
//#include "Tally.h"
//#include "Material.h"
//#include "Isotope.h"
#include "Region.h"
#include "Fissioner.h"
#include "Neutron.h"


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
	void addRegion(Region* region);

	/* Monte Carlo kernel */
	void runMonteCarloSimulation();
//    float computeFuelFuelCollisionProb(int energy_index);
//    float computeModeratorFuelCollisionProb(int energy_index);
};


#endif /* GEOMETRY_H_ */
