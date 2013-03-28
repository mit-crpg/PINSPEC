/*
 * Region.h
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef REGION_H_
#define REGION_H_

#include <vector>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <string.h>
#include "Surface.h"
#include "Material.h"


#define _USE_MATH_DEFINES


typedef enum regionTypes {
	FUEL,
	MODERATOR,
	INFINITE
} regionType;


/**
 * The Region class represents a single dimensional region
 * bounded by two planes, or Surface class objects. The region
 * contains a vector of neutrons which live within it and is filled
 * by a set of Isotope class objects, each with a different number
 * density. The Region class contains all of the physics for moving
 * colliding neutrons
 */

#ifdef __cplusplus
class Region {

private:
	char* _region_name;
	float _volume;
	Material* _material;
	regionType _region_type;

	/* Two region pin cell parameters */
	float _sigma_e;
	float _beta;
	float _alpha1;
	float _alpha2;

	/* Geometry parameters pin cell */
	float _fuel_radius;
	float _pitch;
	float _half_width;

    /* Combinatorial geometry specifications - for HETEROGENEOUS spatial type */
	std::vector<float> _fuel_ring_radii;
	std::vector<float> _moderator_ring_radii;	
    std::vector<Surface*> _bounding_surfaces;
    std::vector<int> _halfspaces;

    bool contains(neutron* neutron);
    bool onBoundary(neutron* neutron);

public:
	Region(char* region_name, regionType type);
	virtual ~Region();
    char* getRegionName();
    float getVolume();
    Material* getMaterial();
    regionType getRegionType();
    bool isFuel();
    bool isModerator();
	bool isInfinite();
	float getFuelRadius();
	float getPitch();

	float getTotalMacroXS(float energy);
	float getTotalMacroXS(int energy_index);
	float getElasticMacroXS(float energy);
	float getElasticMacroXS(int energy_index);
	float getAbsorptionMacroXS(float energy);
	float getAbsorptionMacroXS(int energy_index);
	float getCaptureMacroXS(float energy);
	float getCaptureMacroXS(int energy_index);
	float getFissionMacroXS(float energy);
	float getFissionMacroXS(int energy_index);
	float getTransportMacroXS(float energy);
	float getTransportMacroXS(int energy_index);

    void setMaterial(Material* material);
//    void addBoundingSurface(Surface* surface, int halfspace);

	void setFuelRadius(float radius);
	void setPitch(float pitch);
    void setVolume(float volume);
	void addFuelRingRadius(float radius);
	void addModeratorRingRadius(float radius);

	collisionType collideNeutron(neutron* neut);

};

#endif

#endif /* REGION_H_ */
