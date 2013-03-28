/*
 * Region.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Region.h"


/**
 * Region constructor sets default values
 */
Region::Region(char* region_name, regionType type) {

    _region_name = region_name;
    _region_type = type;

	/* Default two region pin cell parameters */
	_volume = 0.0;
	_fuel_radius = 0.0;
	_pitch = 0.0;
	_half_width = 0.0;

	_sigma_e = 0.0;
	_beta = 0.0;
	_alpha1 = 0.0;
	_alpha2 = 0.0;

}


/**
 * Region destructor does not delete anything since that is left to SWIG
 */
Region::~Region() { }


/**
 * Return the name, if any, of this region as specified by
 * the user
 * @return a character array representing this region's name
 */
char* Region::getRegionName() {
	return _region_name;
}


/**
 * Returns the volume of this Region (cm^3)
 * @return the Region's volume
 */
float Region::getVolume() {
	return _volume;
}


/**
 * Returns a pointer to the Material filling this Region
 * @return a pointer to a Material
 */
Material* Region::getMaterial() {
	return _material;
}


/**
 * Return the type of Region (FUEL, MODERATOR or INFINITE)
 * @return the region type
 */
regionType Region::getRegionType() {
	return _region_type;
}


/**
 * Returns true if this region has been labeled as fuel, false otherwise
 * @return true if fuel, false otherwise
 */
bool Region::isFuel() {
	if (_region_type == FUEL)
		return true;
	else
		return false;
}


/**
 * Returns true if this region has been labeled as moderator, false otherwise
 * @return true if moderator, false otherwise
 */
bool Region::isModerator() {
	if (_region_type == MODERATOR)
		return true;
	else
		return false;
}


/**
 * Returns true if this region has been labeled as infnite, false otherwise
 * @return true if infinite, false otherwise
 */
bool Region::isInfinite() {
	if (_region_type == INFINITE)
		return true;
	else
		return false;
}


/**
 * Returns the radius of the fuel pin if not an INFINITE Region type
 * @return the fuel pin radius
 */
float Region::getFuelRadius() {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot return a fuel pin radius for region %s"
					" which is INFINITE", _region_name);
	return _fuel_radius;
}


/**
 * Returns the pitch of the pin cell if not an INFINITE Region type
 * @return the fuel pin radius
 */
float Region::getPitch() {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot return a pin cell pitch for region %s"
					" which is INFINITE", _region_name);

	return _pitch;
}


float Region::getTotalMacroXS(float energy) {
    return _material->getTotalMacroXS(energy);
}


float Region::getTotalMacroXS(int energy_index) {
    return _material->getTotalMacroXS(energy_index);
}


float Region::getElasticMacroXS(float energy) {
    return _material->getElasticMacroXS(energy);
}


float Region::getElasticMacroXS(int energy_index) {
    return _material->getElasticMacroXS(energy_index);
}


float Region::getAbsorptionMacroXS(float energy) {
    return _material->getAbsorptionMacroXS(energy);
}


float Region::getAbsorptionMacroXS(int energy_index) {
    return _material->getAbsorptionMacroXS(energy_index);
}


float Region::getCaptureMacroXS(float energy) {
    return _material->getCaptureMacroXS(energy);
}

float Region::getCaptureMacroXS(int energy_index) {
    return _material->getCaptureMacroXS(energy_index);
}


float Region::getFissionMacroXS(float energy) {
    return _material->getFissionMacroXS(energy);
}


float Region::getFissionMacroXS(int energy_index) {
    return _material->getFissionMacroXS(energy_index);
}


float Region::getTransportMacroXS(float energy) {
    return _material->getTransportMacroXS(energy);
}


float Region::getTransportMacroXS(int energy_index) {
    return _material->getTransportMacroXS(energy_index);
}


/**
 * Sets this Region's volume (cm^3)
 * @param volume the volume of this Region
 */
void Region::setVolume(float volume) {
	_volume = volume;
}


/**
 * Set the Material filling this Region
 * @param material a pointer to a Material
 */
void Region::setMaterial(Material* material) {
	_material = material;
}


/**
 * Sets the fuel radius for this Region if it is not INFINITE
 * @param radius the fuel pin radius
 */
void Region::setFuelRadius(float radius) {

	_fuel_radius = radius;

    if (_fuel_radius != 0.0) {
        if (_region_type == MODERATOR)
            _volume = _pitch * _pitch - M_PI * _fuel_radius * _fuel_radius;
        else
            _volume = M_PI * _fuel_radius * _fuel_radius;
    }
}


/**
 * Sets the pin cell pitch for this Region if it is not INFINITE
 * @param radius the pin cell pitch
 */
void Region::setPitch(float pitch) {

	_pitch = pitch;
	_half_width = pitch / 2.0;

    if (_fuel_radius != 0.0) {
        if (_region_type == MODERATOR)
            _volume = _pitch * _pitch - M_PI * _fuel_radius * _fuel_radius;
        else
            _volume = M_PI * _fuel_radius * _fuel_radius;
    }
}


/**
 * Adds a new fuel ring radius for this Region if it is not INFINITE
 * @param radius a fuel ring radius
 */
void Region::addFuelRingRadius(float radius) {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot add a fuel pin ring radius for region %s"
					" which is INFINITE", _region_name);
	else if (_region_type == MODERATOR)
		log_printf(ERROR, "Cannot add a fuel ring radius for region %s"
					" which is a MODERATOR type region", _region_name);

	_fuel_ring_radii.push_back(radius);
}



/**
 * Adds a new moderator ring radius for this Region if it is not INFINITE
 * @param radius the fuel pin radius
 */
void Region::addModeratorRingRadius(float radius) {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot add a moderator ring radius for region %s"
					" which is INFINITE", _region_name);
	else if (_region_type == FUEL)
		log_printf(ERROR, "Cannot add a moderator ring radius for region %s"
					" which is a FUEL type region", _region_name);

	_moderator_ring_radii.push_back(radius);
}



/**
 * This method collides a neutron within some Region. It
 * tallies the neutron in any Tallies for this Region and 
 * uses the Region's Material to compute the neutron's 
 * next energy and collision type.
 */
collisionType Region::collideNeutron(neutron* neut) {

	/* Collide the neutron in the Region's Material */
    collisionType type = _material->collideNeutron(neut);

	return type;
}


/**
 * Check if this region contains a certain 2D coordinate location
 * @param x the x-coordinate of the position
 * @param y the y-coordinate of the position
 * @return if contained (true), otherwise (false)
 */
bool Region::contains(neutron* neutron) {

    float x = neutron->_x;
    float y = neutron->_y;

	if  (_region_type == INFINITE)
		return true;
	else {
		float r = pow(x, 2.0) + pow(y, 2.0);
		if (_region_type == FUEL && r < _fuel_radius)
			return true;
		else if (_region_type == MODERATOR && 
					(fabs(x) < _half_width && fabs(y) < _half_width))
			return true;
		else
			return false;			
	}
}


/**
 * Checks if a 2D coordinate location is on a boundary of this Region 
 * @param x the x-coordinate of the position
 * @param y the y-coordinate of the position
 * @return if on boundary (true), otherwise false
 */
bool Region::onBoundary(neutron* neutron) {

    float x = neutron->_x;
    float y = neutron->_y;

	if (_region_type == INFINITE) {
		log_printf(WARNING, "Unable to compute onBoundary method"
					" for region %s since it is INIFINITE", _region_name); 
		return false;
	}
	else {
		float r = pow(x, 2.0) + pow(y, 2.0);
		if (fabs(r - _fuel_radius) < 1E-5)
			return true;
		else if (fabs(x - _half_width) < 1E-5)
			return true;
		else if (fabs(y - _half_width) < 1E-5)
			return true;
		else
			return false;
	}
}

