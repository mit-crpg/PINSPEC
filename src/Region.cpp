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
Region::Region() {

	_region_type = INFINITE;
	_volume = 0.0;

	/* By default the Region is infinite (unbounded) */
	_region_type = INFINITE;
	_spatial_type = HOMOGENEOUS;

	/* Default two region pin cell parameters */
	_sigma_e = 0.0;
	_beta = 0.0;
	_alpha1 = 0.0;
	_alpha2 = 0.0;

	_fuel_radius = 0.0;
	_pitch = 0.0;
	_half_width = 0.0;
}


/**
 * Region destructor
 */
Region::~Region() { }


/**
 * Return the name, if any, of this region as specified by
 * the user
 * @return a character array representing this region's name
 */
char* Region::getRegionName() {
	return _name;
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
 * Return the spatial type of Region (HOMOGENEOUS or HETEROGENEOUS)
 * @return the spatial type
 */
spatialType Region::getSpatialType() {
	return _spatial_type;
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
					" which is INFINITE", _name);
	return _fuel_radius;
}


/**
 * Returns the pitch of the pin cell if not an INFINITE Region type
 * @return the fuel pin radius
 */
float Region::getPitch() {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot return a pin cell pitch for region %s"
					" which is INFINITE", _name);

	return _pitch;
}

/**
 * Sets this Region's name as specified by a user
 * @param region_name the name of this Region
 */
void Region::setRegionName(char* region_name) {
	_name = region_name;
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
 * Set this Region's type (FUEL, MODERATOR, or INFINITE)
 * @param region_type the type of region
 */
void Region::setRegionType(regionType region_type) {
	_region_type = region_type;
}



/**
 * Set this Region's spatial type (HOMOGENEOUS or HETEROGENEOUS)
 * @param spatial_type the spatial type of region
 */
void Region::setSpatialType(spatialType spatial_type) {
	_spatial_type = spatial_type;
}


/**
 * Adds a new Tallies to this region for tallying
 * @param bins a pointer to a Tally class object
 */
void Region::addTally(Tally* tally) {
	_tallies.push_back(tally);
}


/**
 * Sets the escape cross-section, beta, alpha1 and alpha2 parameters
 * used for the two region pin cell simulation
 * @param sigma_e the escape cross-section
 * @param beta Carlvik's beta parameter
 * @param alpha1 Carlvik's alpha1 parameter
 * @param alpha2 Carlvik's alpha2 parameter
 */
void Region::setTwoRegionPinCellParams(float sigma_e, float beta,
											float alpha1, float alpha2) {
	_sigma_e = sigma_e;
	_beta = beta;
	_alpha1 = alpha1;
	_alpha2 = alpha2;
}


/**
 * Sets the fuel radius for this Region if it is not INFINITE
 * @param radius the fuel pin radius
 */
void Region::setFuelRadius(float radius) {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot set a fuel pin radius for region %s"
					" which is INFINITE", _name);

	_fuel_radius = radius;
}


/**
 * Sets the pin cell pitch for this Region if it is not INFINITE
 * @param radius the pin cell pitch
 */
void Region::setPitch(float pitch) {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot set a pin cell pitch for region %s"
					" which is INFINITE", _name);

	_fuel_radius = pitch;
	_half_width = pitch / 2.0;
}


/**
 * Adds a new fuel ring radius for this Region if it is not INFINITE
 * @param radius a fuel ring radius
 */
void Region::addFuelRingRadius(float radius) {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot add a fuel pin ring radius for region %s"
					" which is INFINITE", _name);
	else if (_region_type == MODERATOR)
		log_printf(ERROR, "Cannot add a fuel ring radius for region %s"
					" which is a MODERATOR type region", _name);

	_fuel_ring_radii.push_back(radius);
}



/**
 * Adds a new moderator ring radius for this Region if it is not INFINITE
 * @param radius the fuel pin radius
 */
void Region::addModeratorRingRadius(float radius) {

	if (_region_type == INFINITE)
		log_printf(ERROR, "Cannot add a moderator ring radius for region %s"
					" which is INFINITE", _name);
	else if (_region_type == FUEL)
		log_printf(ERROR, "Cannot add a moderator ring radius for region %s"
					" which is a FUEL type region", _name);

	_moderator_ring_radii.push_back(radius);
}


/**
 * This function computes the two-region fuel-to-fuel collision probability for
 * a two-region pin cell simulation. It uses Carlvik's two-term rational model
 * and assumes that the escape cross-section (_sigma_e), _beta, _alpha1 and
 * _alpha2 have all been set or else it throws an error
 * @param energy_index index into the material's energy grid
 * @return the fuel-to-fuel collision probability at that energy
 */
float Region::computeFuelFuelCollisionProb(int energy_index) {

	float p_ff;

	/* If this is the fuel region, we can compute p_ff directly */
	if (_region_type == FUEL) {

		/* Check that all necessary parameters to compute p_ff have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to compute a fuel-fuel collision "
					"probability for region %s since beta, sigma_e, "
					"alpha1, or alpha2 for this region has not yet been set",
					_name);

		float sigma_tot_fuel = _material->getTotalMacroXS(energy_index);

		p_ff = ((_beta*sigma_tot_fuel) / (_alpha1*_sigma_e + sigma_tot_fuel)) +
			((1.0 - _beta)*sigma_tot_fuel / (_alpha2*_sigma_e + sigma_tot_fuel));
	}

	/* If this is the moderator region, we cannot compute p_ff*/
	else {
		log_printf(ERROR, "Unable to compute fuel-fuel collision "
					"probability for region %s since it is not a "
					"FUEL region", _name);
	}

	return p_ff;
}


/**
 * This function computes the two-region moderator-to-fuel collision
 * probability for a two-region pin cell simulation. It uses Carlvik's
 * two-term rational model and assumes that the escape cross-section
 * (_sigma_e), _beta, _alpha1 and _alpha2 have all been set or else it
 * throws an error
 * @param energy_index index into the material's energy grid
 * @return the moderator-to-fuel collision probability at that energy
 */
float Region::computeModeratorFuelCollisionProb(int energy_index) {

	float p_mf;

	/* If this is the fuel region, we can compute p_mf directly */
	if (_region_type == FUEL) {

		/* Check that all necessary parameters to compute p_mf have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to compute a moderator-fuel collision "
					"probability for region %s since beta, sigma_e, alpha1, "
					"or alpha2 for this region has not yet been set",
					_name);

		float p_ff = computeFuelFuelCollisionProb(energy_index);
		float p_fm = 1.0 - p_ff;

		float tot_sigma_f = _material->getTotalMacroXS(energy_index);
		//NOTE: Need to fix this
//		float tot_sigma_mod = _other_region->getMaterial()->getTotalMacroXS(energy_index);
		float v_mod = _pitch * _pitch - M_PI * _fuel_radius * _fuel_radius;

//		p_mf = p_fm*(tot_sigma_f*_volume) / (tot_sigma_mod*v_mod);
	}

	/* If this is the moderator region, we ask fuel region to compute p_mf */
	else {
		log_printf(ERROR, "Unable to compute moderator-fuel collision "
					"probability for region %s since it is not a FUEL"
					" region", _name);
	}

	return p_mf;
}


/**
 * Clear this region's vector of Tally class object pointers
 */
void Region::clearTallies() {
	_tallies.clear();
}


/**
 * Check if this region contains a certain 2D coordinate location
 * @param x the x-coordinate of the position
 * @param y the y-coordinate of the position
 * @return if contained (true), otherwise (false)
 */
bool Region::contains(float x, float y) {

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
bool Region::onBoundary(float x, float y) {

	if (_region_type == INFINITE) {
		log_printf(WARNING, "Unable to compute onBoundary method"
					" for region %s since it is INIFINITE", _name); 
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
