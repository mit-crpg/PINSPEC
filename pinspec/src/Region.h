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
    /** A fuel region */
    FUEL,
    /** A moderator region */
    MODERATOR,
    /** An infinite medium region */
    INFINITE
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

private:
    /** The region's name */
    char* _region_name;
    /** The volume occupied by the region in 2D space */
    float _volume;
    /** A pointer to the material filling the region */
    Material* _material;
    /** The type of region (INFINITE, MODERATOR, FUEL, etc.) */
    regionType _region_type;

  /* Geometry parameters pin cell */
    /** The radius of the fuel */
    float _fuel_radius;
    /** The pin cell pitch */
    float _pitch;
    /** Half of the pin cell pitch */
    float _half_width;
    /** The squared geometric buckling */
    float _buckling_squared;

    bool contains(neutron* neutron);
    bool onBoundary(neutron* neutron);

public:
    Region(char* region_name, regionType type);
    virtual ~Region();
    char* getRegionName();
    float getVolume();
    Material* getMaterial();
    bool containsIsotope(Isotope* isotope);
    regionType getRegionType();
    bool isFuel();
    bool isModerator();
    bool isInfinite();
    float getFuelRadius();
    float getPitch();
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
    void setFuelRadius(float radius);
    void setPitch(float pitch);
    void setVolume(float volume);
    void setBucklingSquared(float buckling_squared);

    void collideNeutron(neutron* neutron);
};


#endif /* REGION_H_ */
