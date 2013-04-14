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
    InfiniteMediumRegion* _infinite_medium;
    /** FUEL type region if the geometry is HOMOGENEOUS_EQUIVALANCE */
    EquivalenceRegion* _fuel;
    /** MODERATOR type region if the geometry is HOMOGENEOUS_EQUIVALANCE */
    EquivalenceRegion* _moderator;
    /** A container of BOUNDED type regions if the geometry is HETEROGENEOUS */
    std::vector<BoundedRegion*> _regions;

    /** The fuel pin radius for a heterogeneous-homogeneous equivalent
     *  geometry */
    float _fuel_radius;
    /** The pin cell pitch for a heterogenous-homogeneous equivalent
     *  geometry */
    float _pitch;
    /** The square of the geometric buckling */
    float _buckling_squared;

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

    /** The number of first flight collision probabilities */
    int _num_prob;
    /** The energies on a uniform lethargy grid for which the first flight 
     * collision probabilities are defined */
    float* _prob_energies;
    /** The first flight fuel-to-fuel collision probabilities */
    float* _prob_ff;
    /** The first flight moderator-to-fuel collision probabilities */
    float* _prob_mf;

    /** The fissioner used to sample new neutron fission emission energies */
    Fissioner* _fissioner;
    /** A 3D spherical radius within which to sample random source sites */
    float _source_sampling_radius;

    void initializeProbModFuelRatios();
    void findContainingRegion(neutron* neutron);
    bool contains(neutron* neutron);

public:
    Geometry(spatialType spatial_type);
    virtual ~Geometry();
	
    int getNumNeutronsPerBatch();
    int getTotalNumNeutrons();
    int getNumBatches();
    int getNumThreads();
    spatialType getSpatialType();
    float getBucklingSquared();
    float getVolume();
    float getSourceSamplingRadius();

    void setSourceSamplingRadius(float radius);			
    void setNeutronsPerBatch(int num_neutrons_per_batch);
    void setNumBatches(int num_batches);
    void setNumThreads(int num_threads);
    void setSpatialType(spatialType spatial_type);
    void setFuelPinRadius(float radius);
    void setPinCellPitch(float pitch);
    void setDancoffFactor(float dancoff);
    void addRegion(Region* region);
    void setBucklingSquared(float buckling_squared);

    void runMonteCarloSimulation();
    void initializeSourceNeutron(neutron* neutron);
};


#endif /* GEOMETRY_H_ */
