/**
 * @file Region.h
 * @brief The Region class.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 15, 2013
 */

#ifndef REGION_H_
#define REGION_H_

#ifdef __cplusplus
#include <vector>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <string.h>
#include <utility>
#include "Surface.h"
#include "Material.h"
#endif

#define _USE_MATH_DEFINES

/**
 * @enum regionTypes
 * @brief Spatial Types of regions
 */

/**
 * @var regionType
 * @brief Region spatial type
 */
typedef enum regionTypes {
    /** An infinite medium region */
    INFINITE_MEDIUM,
    /** An equivalence theory-based fuel region */
    EQUIVALENT_FUEL,
    /** An equivalence theory-based fuel moderator */
    EQUIVALENT_MODERATOR,
    /** A fully heterogeneous circular fuel region bounded by surfaces */
    BOUNDED_FUEL,
    /** A fully heterogeneous pin cell moderator region bounded by surfaces */
    BOUNDED_MODERATOR,
    /** A fully heterogenous generalized region bounded by surfaces */
    BOUNDED_GENERAL
} regionType;



/**
 * @class Region Region.h "src/pinspec/Region.h"
 * @brief The region class represents a region in 2D space.
 * @details The region class is a superclass allowing for subclasses 
 *          representing infinited media, heterogeneous/homogeneous
 *          equivalance fuel/moderator regions, or even heterogeneous
 *          regions bounded by 2D quadratic surfaces, filled by a material.
 */
class Region {

protected:
    /** The region's name */
    char* _region_name;
    /** A static class variable to generate a UID for each new region */
    static int _n;
    /** The region's unique identifier */
    int _uid;
    /** A pointer to the material filling the region */
    Material* _material;
    /** The type of region (INFINITE, EQUIVALNENCE, or BOUNDED) */
    regionType _region_type;
    /** The squared geometric buckling */
    float _buckling_squared;
    /** The volume occupied by the region in 2D space */
    float _volume;
    /** The random number seed */
    unsigned int _seed;

public:
    Region(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    virtual ~Region();
    char* getName();
    int getUid() const;
    Material* getMaterial();
    bool containsIsotope(Isotope* isotope);
    regionType getRegionType();
    float getVolume();
    float getBucklingSquared();

    float getTotalMacroXS(float energy);
    float getTotalMacroXS(int energy_index);
    float getTotalMicroXS(float energy);
    float getTotalMicroXS(int energy_index);
	
    float getElasticMacroXS(float energy);
    float getElasticMacroXS(int energy_index);
    float getElasticMicroXS(float energy);
    float getElasticMicroXS(int energy_index);
	
    float getAbsorptionMacroXS(float energy);
    float getAbsorptionMacroXS(int energy_index);
    float getAbsorptionMicroXS(float energy);
    float getAbsorptionMicroXS(int energy_index);
	
    float getCaptureMacroXS(float energy);
    float getCaptureMacroXS(int energy_index);
    float getCaptureMicroXS(float energy);
    float getCaptureMicroXS(int energy_index);
	
    float getFissionMacroXS(float energy);
    float getFissionMacroXS(int energy_index);
    float getFissionMicroXS(float energy);
    float getFissionMicroXS(int energy_index);
	
    float getTransportMicroXS(float energy);
    float getTransportMicroXS(int energy_index);
    float getTransportMacroXS(float energy);
    float getTransportMacroXS(int energy_index);

    void setMaterial(Material* material);
    void setVolume(float volume);
    void setBucklingSquared(float buckling_squared);
    void setRandomNumberSeed(unsigned int seed);
    void initializeRandomNumberGenerator();

    /**
     * @brief This method collides a neutron within the region.
     * @details This method encapsulates all of the neutron scattering physics
     *          which is further encapsulated by the material and isotope 
     *          classes. This method must be implemented for each region
     *          subclass type.
     * @param neutron the neutron of interest
     */
    virtual void collideNeutron(neutron* neutron) =0;
};


/**
 * @class InfiniteMediumRegion Region.h "pinspec/src/Region.h"
 * @brief The InfiniteMediumRegion is a homogenized region treated without
 *        geometric effects.
 * @details Only a single InfiniteMedium Region should be defined for a given
 *          PINSPEC simulation. An infinite medium spectral calculation in
 *          in PINSPEC neglects the location of particles and only treats
 *          neutron elastic and thermal capture, scattering and fission events.
 */
class InfiniteMediumRegion: public Region {

public:
    InfiniteMediumRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    virtual ~InfiniteMediumRegion() { };

    void collideNeutron(neutron* neutron);
};


/**
 * @class EquivalenceRegion Region.h "pinspec/src/Region.h"
 * @brief The Equivalence Region is a fuel or moderator region treated using
 *        the heterogenous-homogeneous equivalence theory. 
 * @details In particular, this class uses Carlvik's Rational Approximation 
 *        for first flight collision probabilities. This heterogeneous -
 *        homogeneous approximation allows this class to treat a heterogeneous 
 *        pin cell simply as a homogeneous fuel and homogeneous moderator 
 *        regions with probabilities for a particle to travel from one region 
 *        to another at a given energy. This avoids the need for explicit ray 
 *        tracing of particle trajectories across the geometry.
 */
class EquivalenceRegion: public Region {

private:
    /** The radius of the fuel */
    float _fuel_radius;
    /** The pin cell pitch */
    float _pitch;
    /** Half of the pin cell pitch */
    float _half_width;

    /** The other region (ie, if this region is an EQUIVALENT_FUEL type region,
     *the other region is EQUIVALENT_MODERATOR) */
    EquivalenceRegion* _other_region;

    /** First flight fuel-to-fuel collision probabilities */
    float* _prob_ff;
    /** First flight fuel-to-fuel collision probabilities */
    float* _prob_mf;
    /** Energies (eV) for first flight collision probabilities */
    float* _prob_energies;
    /** The number of first flight collision probabilities */
    int _num_prob;
    /** Lowest lethargy for the first flight collision probabilities */
    float _start_lethargy;
    /** Highest lethargy for the first flight collision probabilities */
    float _end_lethargy;
    /** Space between lethargy for the first flight collision probabilities*/
    float _delta_lethargy;

public:
    EquivalenceRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    virtual ~EquivalenceRegion() { };
    float getFuelPinRadius();
    float getPinCellPitch();
    int getEnergyGridIndex(float lethargy);
    bool isFuel();
    bool isModerator();

    void setFirstFlightCollProb(float* prob_ff, float* prob_mf, 
				float* prob_energies, int num_prob);
    void setOtherRegion(EquivalenceRegion* region);
    void setFuelPinRadius(float radius);
    void setPinCellPitch(float pitch);

    float computeFuelFuelCollsionProb(neutron* neutron);
    float computeModeratorFuelCollisionProb(neutron* neutron);
    void collideNeutron(neutron* neutron);
};


/**
 * @class EquivalenceFuelRegion Region.h "pinspec/src/Region.h"
 * @brief The EquivalenceFuelRegion is a fuel region treated using the
 *        heterogenous-homogeneous equivalence theory. 
 * @details In particular, this class uses Carlvik's Rational Approximation 
 *        for first flight collision probabilities. This heterogeneous -
 *        homogeneous approximation allows this class to treat a heterogeneous 
 *        pin cell simply as a homogeneous fuel and homogeneous moderator 
 *        regions with probabilities for a particle to travel from one region 
 *        to another at a given energy. This avoids the need for explicit ray 
 *        tracing of particle trajectories across the geometry.
 */
class EquivalenceFuelRegion: public EquivalenceRegion {

public:
    EquivalenceFuelRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    ~EquivalenceFuelRegion() { };
};


/**
 * @class EquivalenceModeratorRegion Region.h "pinspec/src/Region.h"
 * @brief The EquivalenceModeratorRegion is a fuel region treated using the
 *        heterogenous-homogeneous equivalence theory. 
 * @details In particular, this class uses Carlvik's Rational Approximation 
 *        for first flight collision probabilities. This heterogeneous -
 *        homogeneous approximation allows this class to treat a heterogeneous 
 *        pin cell simply as a homogeneous fuel and homogeneous moderator 
 *        regions with probabilities for a particle to travel from one region 
 *        to another at a given energy. This avoids the need for explicit ray 
 *        tracing of particle trajectories across the geometry.
 */
class EquivalenceModeratorRegion: public EquivalenceRegion {

public:
    EquivalenceModeratorRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    ~EquivalenceModeratorRegion() { };
};


/**
 * @class BoundedRegion Region.h "pinspec/src/Region.h"
 * @brief The BoundedRegion is a region bounded by surfaces using a constructive
 *        solid geometry formulation.
 * @details A BoundedRegion is a space that is contained by the intersection 
 *          of halfspaces of surface boundaries using a constructive solid 
 *          geometry formulation. It is up to the user to ensure that the
 *          surfaces added to a BoundedRegion do indeed form a closed region.
 */
class BoundedRegion: public Region {

private:
    /** A container of pointers to the surfaces bounding this region */
    std::vector< std::pair<int, Surface*> > _surfaces;

public:
    BoundedRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    virtual ~BoundedRegion() { };
    void addBoundingSurface(int halfspace, Surface* surface);
    void removeBoundingSurface(int halfspace, Surface* surface);

    bool contains(neutron* neutron);
    bool contains(float x, float y, float z);
    bool onBoundary(neutron* neutron);
    float computeParametrizedDistance(neutron* neutron);

    void collideNeutron(neutron* neutron);
};




/**
 * @class BoundedFuelRegion Region.h "pinspec/src/Region.h"
 * @brief The BoundedFuelRegion is a circular fuel pin region bounded by 
 *        the interior halfspace of a circle surface.
 * @details A BoundedFuelRegion is a specialized bounded surface object
 *          which allows the user to simply and easily subdivide the
 *          region into equal areaa rings. Each ring will be a cloned version 
 *          of the original region but will contain different bounding
 *          surfaces and new, cloned tallies. This allows for a more 
 *          descriptive picture of the flux gradient across a fuel pin.
 */
class BoundedFuelRegion: public BoundedRegion {

public:
    BoundedFuelRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    ~BoundedFuelRegion() { };

    void ringify(int num_rings);
};


/**
 * @class BoundedModeratorRegion Region.h "pinspec/src/Region.h"
 * @brief The BoundedModeratorRegion is a the space defined four surfaces
 *        bounding a square pin minus a circular fuel pin region.
 * @details A BoundedModeratorRegion is a specialized bounded surface object
 *          which allows the user to simply and easily subdivide the
 *          region into equal area rings. Each ring will be a cloned version 
 *          of the original region but will contain different bounding
 *          surfaces and new, cloned tallies. This allows for a more 
 *          descriptive picture of the flux gradient within the moderator.
 */
class BoundedModeratorRegion: public BoundedRegion {

public:
    BoundedModeratorRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    ~BoundedModeratorRegion() { };

    void ringify(int num_rings);
};


/**
 * @class BoundedGeneralRegion Region.h "pinspec/src/Region.h"
 * @brief The BoundedGeneralRegion is a the space defined by the intersection
 *        of some number surface halfspaces.
 * @details This type of region is useful for a region that is more 
 *          complicated than the basic fuel or moderator pin cell, and
 *          allows the user to apply their own discretion to create a
 *          the bounded region.
 */
class BoundedGeneralRegion: public BoundedRegion {

public:
    BoundedGeneralRegion(const char* region_name=(char*)"");
    /**
     * @brief Empty destructor allows SWIG to cleanup memory for surfaces.
     */
    ~BoundedGeneralRegion() { };
};

#endif /* REGION_H_ */

