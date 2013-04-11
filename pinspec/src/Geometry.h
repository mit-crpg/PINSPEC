/**
 * @file Geometry.h
 * @brief The Geometry class.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 6, 2013.
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
#endif


/**
 * @enum spatialTypes
 * @brief The spatial types for the geometry
 */

/**
 * @var spatialType
 * @brief A spatial type for the geometry
 */
typedef enum spatialTypes {
    /** An infinite homogeneous medium */
    INFINITE_HOMOGENEOUS,
    /** A heterogeneous-homogeneous equivalent geometry */
    HOMOGENEOUS_EQUIVALENCE,
    /** A heterogenous geometry built of regions and surfaces */
    HETEROGENEOUS,
} spatialType;


/**
 * @class Geometry Geometry.h "pinspec/src/Geometry.h"
 * @brief The Geometry represents the highest level entity in which a neutron
 *        may reside during a PINSPEC simulaiton. The geometry consists of
 *        one or more regions and controls the highest level Monte Carlo
 *        kernel.
 */
class Geometry {

private:

    /** The number of neutrons per batch */
    int _num_neutrons_per_batch;
    /** The number of batches */
    int _num_batches;
    /** The number of threads */
    int _num_threads;

    /** The spatial type for the geometry */
    spatialType _spatial_type;
    /** INFINITE type region if the geometery is INFINITE_HOMOGENEOUS */
    Region* _infinite_medium;
    /** FUEL type region if the geometry is HETEROGENEOUS or 
     *  HOMOGENEOUS_EQUIVALANCE */
    Region* _fuel;
    /** MODERATOR type region if the geometry is HETEROGENEOUS or 
     *  HOMOGENEOUS_EQUIVALANCE */
    Region* _moderator;
    /** An array of neutrons */
    neutron* _neutrons;
    /** The fissioner used to sample new neutron fission emission energies */
    Fissioner* _fissioner;

    /** The user-specified dancoff factor */    
    float _dancoff;
    /** The user-specified escape cross-section */
    float _sigma_e;
    /** The user-specified beta value for Carlvik's rational approximation */
    float _beta;
    /** The user-specified alpha1 value for Carlvik's rational approximation */
    float _alpha1;
    /** The user-specified alpha2 value for Carlvik's rational approximation */
    float _alpha2;

    /** The square of the geometric buckling */
    float _buckling_squared;

    /** The number of moderator to fuel cross-section ratios */
    int _num_ratios;
    /** An array of the moderator to fuel cross-section ratios */
    float* _pmf_ratios;
    /** The spacing type between bins (EQUAL, LOGARITHMIC, OTHER) */
    binSpacingTypes _scale_type;
    /** Lowest energy for the moderator to fuel cross-section ratios */
    float _start_energy;
    /** Highest energy for the moderator to fuel cross-section ratios */
    float _end_energy;
    /** Space between energies for the moderator to fuel cross-section ratios */
    float _delta_energy;

    void initializeProbModFuelRatios();
    int getEnergyGridIndex(float energy) const;
    float computeFuelFuelCollisionProb(neutron* neutron);
    float computeModeratorFuelCollisionProb(neutron* neutron);

public:
    Geometry();
    virtual ~Geometry();
	
    int getNumNeutronsPerBatch();
    int getTotalNumNeutrons();
    int getNumBatches();
    int getNumThreads();
    spatialType getSpatialType();
    float getBucklingSquared();
    float getVolume();
			
    void setNeutronsPerBatch(int num_neutrons_per_batch);
    void setNumBatches(int num_batches);
    void setNumThreads(int num_threads);
    void setSpatialType(spatialType spatial_type);
    void setDancoffFactor(float dancoff);
    void addRegion(Region* region);
    void setBucklingSquared(float buckling_squared);

    void runMonteCarloSimulation();
};


/**
 * @brief This method returns the index for a certain energy (eV) into
 *        the uniform lethargy grid of the geometry's first flight
 *        probability of travel from moderator to fuel \f$ p_{mf} \f$ ratios
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


#endif /* GEOMETRY_H_ */
