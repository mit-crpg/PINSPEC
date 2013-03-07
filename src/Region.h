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
#include <omp.h>
#include "Material.h"
#include "Isotope.h"
#include "Neutron.h"
#include "Tally.h"

#define _USE_MATH_DEFINES


typedef enum regionTypes {
	FUEL,
	MODERATOR,
	INFINITE
} regionType;


typedef enum spatialTypes {
	NONE,
	HETEROGENEOUS,
	HOMOGENEOUS
} spatialType;


	/**
 * The Region class represents a single dimensional region
 * bounded by two planes, or Surface class objects. The region
 * contains a vector of neutrons which live within it and is filled
 * by a set of Isotope class objects, each with a different number
 * density. The Region class contains all of the physics for moving
 * colliding neutrons
 */
class Region {

private:
	char* _name;
	float _volume;
	Material* _material;
	regionType _region_type;
	spatialType _spatial_type;

	/* Tallies */
	std::vector<Tally*> _tallies;

	/* Two region pin cell parameters */
	float _sigma_e;
	float _beta;
	float _alpha1;
	float _alpha2;

	/* Geometry parameters pin cell - needed for HETEROGENOUS spatial type */
	float _fuel_radius;
	float _pitch;
	float _half_width;
	std::vector<float> _fuel_ring_radii;
	std::vector<float> _moderator_ring_radii;	

public:
	Region();
	virtual ~Region();
    char* getRegionName();
    float getVolume();
    Material* getMaterial();
    regionType getRegionType();
    spatialType getSpatialType();
    bool isFuel();
    bool isModerator();
	bool isInfinite();
	float getFuelRadius();
	float getPitch();

    void setRegionName(char* region_name);
    void setVolume(float volume);
    void setMaterial(Material* material);
	void setRegionType(regionType region_type);
	void setSpatialType(spatialType spatial_type);
    void addTally(Tally* bins);
    void setTwoRegionPinCellParams(float sigma_e, float beta, 
									float alpha1, float alpha2);
	void setFuelRadius(float radius);
	void setPitch(float pitch);
	void addFuelRingRadius(float radius);
	void addModeratorRingRadius(float radius);

    float computeFuelFuelCollisionProb(int energy_index);
    float computeModeratorFuelCollisionProb(int energy_index);
    void clearTallies();
    bool contains(float x, float y);
    bool onBoundary(float x, float y);
};

#endif /* REGION_H_ */
