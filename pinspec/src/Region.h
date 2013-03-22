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
class Region {

private:
	char* _region_name;
	float _volume;
	Material* _material;
	regionType _region_type;

	/* Geometry parameters pin cell */
	float _fuel_radius;
	float _pitch;
	float _half_width;
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
