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
    _material = NULL;

    /* Default volume */
    if (_region_type == INFINITE)
    	_volume = 1.0;
    else
        _volume = 0.0;

	/* Default two region pin cell parameters */
	_fuel_radius = 0.0;
	_pitch = 0.0;
	_half_width = 0.0;
	_buckling_squared = 0.0;
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


bool Region::containsIsotope(Isotope* isotope) {
	return _material->containsIsotope(isotope);
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


float Region::getBucklingSquared() {
	return _buckling_squared;
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

    if (_material != NULL)
        _material->incrementVolume(_volume);
}


/**
 * Set the Material filling this Region
 * @param material a pointer to a Material
 */
void Region::setMaterial(Material* material) {
	_material = material;
    _material->incrementVolume(_volume);
}


/**
 * Sets the fuel radius for this Region if it is not INFINITE
 * @param radius the fuel pin radius
 */
void Region::setFuelRadius(float radius) {

	_fuel_radius = radius;

    if (_pitch != 0.0) {
        if (_region_type == MODERATOR)
            _volume = _pitch * _pitch - M_PI * _fuel_radius * _fuel_radius;
        else
            _volume = M_PI * _fuel_radius * _fuel_radius;
    }

    if (_material != NULL)
        _material->incrementVolume(_volume);
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

    if (_material != NULL)
        _material->incrementVolume(_volume);
    }
}


void Region::setBucklingSquared(float buckling_squared) {
	_buckling_squared = buckling_squared;
	_material->setBucklingSquared(_buckling_squared);
}


/**
 * This method collides a neutron within some Region. It
 * tallies the neutron in any Tallies for this Region and 
 * uses the Region's Material to compute the neutron's 
 * next energy and collision type.
 */
void Region::collideNeutron(neutron* neutron) {


	if (_material == NULL){
		log_printf(ERROR, "Region %s must have material to"
				" collide neutron", _region_name);

	}

	/* Collide the neutron in the Region's Material */
	neutron->_material = _material;
    _material->collideNeutron(neutron);

	return;
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
		float r = pow((pow(x, 2.0) + pow(y, 2.0)), 0.5);
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
		float r = pow((pow(x, 2.0) + pow(y, 2.0)), 0.5);
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

