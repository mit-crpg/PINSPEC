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
 * Region destructor does not delete anything since that is left to SWIG
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
 * Adds a new Tallies to this region for tallying
 * @param bins a pointer to a Tally class object
 */
void Region::addTally(Tally* tally) {

    if (tally->getTallyDomainType() != REGION)
        log_printf(ERROR, "Unable to add Tally %s to Region %s since the Tally"
                        " is not for an REGION tally domain", 
                                        tally->getTallyName(), _name);
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
 * Sets the number of batches for each of the Tallies inside of this Region
 * @param num_batches the number of batches
 */
void Region::setNumBatches(int num_batches) {

    /* Set the number of batches for each Tally inside of this Region */
    std::vector<Tally*>::iterator iter;
	for (iter = _tallies.begin(); iter != _tallies.end(); iter ++) {
        (*iter)->setNumBatches(num_batches);
    }

    /* Set the number of batches for each of the Tallies inside the Material */
    _material->setNumBatches(num_batches);

    return;
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
 * This method collides a neutron within some Region. It
 * tallies the neutron in any Tallies for this Region and 
 * uses the Region's Material to compute the neutron's 
 * next energy and collision type.
 */
void Region::collideNeutron(neutron* neut) {

	/* If this Region contains any tallies, tally neutron 
	 * before collision
	 */
    int batch_num = neut->_batch_num;
    float sample = neut->_energy;
	float total_xs = _material->getTotalMacroXS(sample);

	/* Collide the neutron in the Region's Material */
    collisionType type = _material->collideNeutron(neut);

    std::vector<Tally*>::iterator iter;
	for (iter = _tallies.begin(); iter != _tallies.end(); iter ++) {
	    Tally *tally = *iter;
	    tallyType tally_type = tally->getTallyType();
	    switch (tally_type) {
	    case FLUX:
			tally->weightedTally(sample, 1.0 / total_xs, batch_num);
	    case COLLISION_RATE:
		    tally->weightedTally(sample, 1.0, batch_num);
	    case ELASTIC_RATE:
		if (type == ELASTIC)
		    tally->weightedTally(sample, 
			_material->getElasticMacroXS(sample) / total_xs, batch_num);
	    case ABSORPTION_RATE:
		if (type == CAPTURE || type == FISSION)
		    tally->weightedTally(sample, 
			_material->getAbsorptionMacroXS(sample) / total_xs, batch_num);
	    case CAPTURE_RATE:
		if (type == CAPTURE)
		    tally->weightedTally(sample, 
			_material->getCaptureMacroXS(sample) / total_xs, batch_num);
	    case FISSION_RATE:
		if (type == FISSION)
		    tally->weightedTally(sample, 
			_material->getFissionMacroXS(sample) / total_xs, batch_num);
	    case TRANSPORT_RATE:
		if (type == ELASTIC)
		    tally->weightedTally(sample, 
			_material->getTransportMacroXS(sample) / total_xs, batch_num);
	    case DIFFUSION_RATE: /* FIXME */
		if (type == ELASTIC)
		    tally->weightedTally(sample, 
			 1.0 / (3.0 * _material->getTransportMacroXS(sample) * total_xs), 
																batch_num); 
	    case LEAKAGE_RATE:; /* FIXME */
	    }
	}

	return;
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


/**
 * Calls each of the Tally class objects in the Region to compute
 * their batch-based statistics from the tallies
 */
void Region::computeBatchStatistics() {

    /* Compute statistics for each of this Region's Tallies */
    std::vector<Tally*>::iterator iter;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter)
        (*iter)->computeBatchStatistics();

    /* Compute statistics for the Material's Tallies */
    _material->computeBatchStatistics();

    return;
}


/**
 * Calls each of the Tally class objects in the Region to compute
 * their batch-based statistics from the tallies
 */
void Region::computeScaledBatchStatistics(float scale_factor) {

    /* Account for the volume of the region if it is not infinite */
	if (_region_type != INFINITE)
        scale_factor *= _volume;

    /* Compute statistics for each of this Region's Tallies */
    std::vector<Tally*>::iterator iter;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter)
        (*iter)->computeScaledBatchStatistics(scale_factor);

    /* Compute statistics for the Material's Tallies */
    _material->computeScaledBatchStatistics(scale_factor);

    return;
}


/**
 * Calls each of the Tally class objects in the Region to output
 * their tallies and statistics to output files.
 * @param directory the directory to write batch statistics files
 * @param suffix a string to attach to the end of each filename
 */
void Region::outputBatchStatistics(char* directory, char* suffix) {

    /* Output statistics for each of this Region's Tallies */
    std::vector<Tally*>::iterator iter;
    std::string filename;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter) {
        filename = std::string(directory) + _name + "_statistics_" 
                                        + suffix + ".txt";
        (*iter)->outputBatchStatistics(filename.c_str());
    }

    /* Output statistics for the Materials Tallies */
    _material->outputBatchStatistics(directory, suffix);

    return;
}

