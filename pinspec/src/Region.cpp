#include "Region.h"


int Region::_n = 1;


/**
 * @brief Region constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_name the (optional) name of the region 
 */
Region::Region(const char* region_name) {

    int length = strlen(region_name);
    _region_name = new char[length];
   
    for (int i=0; i <= length; i++)
      _region_name[i] = region_name[i];

    _uid = _n;
    _n++;
    _material = NULL;
    _volume = 1.0;
    _buckling_squared = 0.0;
}


/**
 * @brief Destructor lets SWIG delete the materials and isotopes during 
 *        garbage collection.
 */
Region::~Region() { 
    delete [] _region_name;
}


/**
 * @brief Return the name of the region.
 * @return a character array representing this region's name
 */
char* Region::getName() {
    return _region_name;
}


/**
 * @brief Returns the unique ID auto-generated for the region.
 * @return a unique ID for the region
 */
int Region::getUid() const {
    return _uid;
}


/**
 * @brief Returns a pointer to the material filling this region.
 * @return a pointer to the material filling the region
 */
Material* Region::getMaterial() {
    return _material;
}


/**
 * @brief Determines whether this region contains a particular isotope.
 * @param isotope the isotope of interest
 * @return true if the region contains the isotope; otherwise false
 */
bool Region::containsIsotope(Isotope* isotope) {
    return _material->containsIsotope(isotope);
}


/**
 * @brief Return the type of region.
 * @return the region type (INFINITE, EQUIVALENT_FUEL, etc.)
 */
regionType Region::getRegionType() {
    return _region_type;
}


/**
 * @brief returns the volume of this Region \f$ (cm^3) \f$.
 * @return the region's volume
 */
float Region::getVolume() {
    return _volume;
}


/**
 * @brief Returns the squared geometric buckling.
 * @return the geometric buckling squared
 */
float Region::getBucklingSquared() {
    return _buckling_squared;
}


/**
 * @brief Computes and returns the total macroscopic cross-section in the region
 *        at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total macroscopic cross-section \f$ (cm^{-1}) \f$
 */ 
float Region::getTotalMacroXS(float energy) {
    return _material->getTotalMacroXS(energy);
}


/**
 * @brief Computes and returns the total macroscopic cross-section in the region
 *        at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the total macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Region::getTotalMacroXS(int energy_index) {
    return _material->getTotalMacroXS(energy_index);
}


/**
 * @brief Computes and returns the total microscopic cross-section in the region
 *        at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total microscopic cross-section
 */ 
float Region::getTotalMicroXS(float energy) {
    return _material->getTotalMicroXS(energy);
}


/**
 * @brief Computes and returns the total microscopic cross-section in the region
 *        at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the total microscopic cross-section
 */
float Region::getTotalMicroXS(int energy_index) {
    return _material->getTotalMicroXS(energy_index);
}


/**
 * @brief Computes and returns the macroscopic elastic scattering 
 *        cross-section in the region at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the macroscopic elastic scattering cross-section \f$ (cm^{-1}) \f$
 */
float Region::getElasticMacroXS(float energy) {
    return _material->getElasticMacroXS(energy);
}


/**
 * @brief Computes and returns the macroscopic elastic scattering cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the total macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Region::getElasticMacroXS(int energy_index) {
    return _material->getElasticMacroXS(energy_index);
}


/**
 * @brief Computes and returns the microscopic elastic scattering  
 *        cross-section in the region at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the microscopic elastic scattering cross-section
 */ 
float Region::getElasticMicroXS(float energy) {
    return _material->getElasticMicroXS(energy);
}  


/**
 * @brief Computes and returns the microscopic elastic scattering cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the microscopic elastic scattering cross-section
 */
float Region::getElasticMicroXS(int energy_index) {
    return _material->getElasticMicroXS(energy_index);
}


/**
 * @brief Computes and returns the macroscopic absorption cross-section
 *        in the region at energy (eV).
 * @param energy the energy of interest (eV)
 * @return the macroscopic absorpotion cross-section \f$ (cm^{-1}) \f$
 */
float Region::getAbsorptionMacroXS(float energy) {
    return _material->getAbsorptionMacroXS(energy);
}


/**
 * @brief Computes and returns the macroscopic absorption cross-section 
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the macroscopic absorption cross-section \f$ (cm^{-1}) \f$
 */
float Region::getAbsorptionMacroXS(int energy_index) {
    return _material->getAbsorptionMacroXS(energy_index);
}


/**
 * @brief Computes and returns the microscopic absorption cross-section
 *        in the region at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the microscopic absorption cross-section
 */
float Region::getAbsorptionMicroXS(float energy) {
    return _material->getAbsorptionMicroXS(energy);
}


/**
 * @brief Computes and returns the microscopic absorption cross-section 
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the microscopic absorption cross-section
 */
float Region::getAbsorptionMicroXS(int energy_index) {
    return _material->getAbsorptionMicroXS(energy_index);
}


/**
 * @brief Computes and returns the macroscopic capture cross-section 
 *        in the region at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the macroscopic capture cross-section \f$ (cm^{-1}) \f$
 */
float Region::getCaptureMacroXS(float energy) {
    return _material->getCaptureMacroXS(energy);
}


/**
 * @brief Computes and returns the macroscopic capture cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the macroscopic capture cross-section \f$ (cm^{-1}) \f$
 */
float Region::getCaptureMacroXS(int energy_index) {
    return _material->getCaptureMacroXS(energy_index);
}


/**
 * @brief Computes and returns the microscopic capture cross-section
 *        in the region at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the microscopic capture cross-section
 */
float Region::getCaptureMicroXS(float energy) {
    return _material->getCaptureMicroXS(energy);
}


/**
 * @brief Computes and returns the microscopic capture cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the microscopic capture cross-section
 */
float Region::getCaptureMicroXS(int energy_index) {
    return _material->getCaptureMicroXS(energy_index);
}


/**
 * @brief Computes and returns the macroscopic fission cross-section 
 *        in the region at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the macroscopic fission cross-section \f$ (cm^{-1}) \f$
 */
float Region::getFissionMacroXS(float energy) {
    return _material->getFissionMacroXS(energy);
}


/**
 * @brief Computes and returns the macroscopic fission cross-section 
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the macroscopic capture cross-section \f$ (cm^{-1}) \f$
 */
float Region::getFissionMacroXS(int energy_index) {
    return _material->getFissionMacroXS(energy_index);
}


/**
 * @brief Computes and returns the microscopic fission cross-section
 *        in the region at some energy (eV)
 * @param energy the energy of interest (eV)
 * @return the microscopic fission cross-section
 */
float Region::getFissionMicroXS(float energy) {
    return _material->getFissionMicroXS(energy);
}


/**
 * @brief Computes and returns the microscopic fission cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the microscopic fission cross-section
 */
float Region::getFissionMicroXS(int energy_index) {
    return _material->getFissionMicroXS(energy_index);
}


/**
 * @brief Computes and returns the microscopic transport cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy the energy of interest
 * @return the microscopic transport cross-section
 */
float Region::getTransportMicroXS(float energy) {
    return _material->getTransportMicroXS(energy);
}


/**
 * @brief Computes and returns the microscopic transport cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the microscopic transport cross-section
 */
float Region::getTransportMicroXS(int energy_index) {
    return _material->getTransportMicroXS(energy_index);
}


/**
 * @brief Computes and returns the macroscopic transport cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy the energy of interest (eV)
 * @return the macroscopic transport cross-section
 */
float Region::getTransportMacroXS(float energy) {
    return _material->getTransportMacroXS(energy);
}


/**
 * @brief Computes and returns the macroscopic transport cross-section
 *        in the region at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the macroscopic transport cross-section
 */
float Region::getTransportMacroXS(int energy_index) {
    return _material->getTransportMacroXS(energy_index);
}


/**
 * @brief Sets the volume for this region \f$ (cm^3) \f$.
 * @param volume the volume occuppied by this region
 */
void Region::setVolume(float volume) {
    _volume = volume;

    if (_material != NULL)
        _material->incrementVolume(_volume);
}


/**
 * @brief Set the material filling this region.
 * @details This method also increments the volume for the material by the
 *          volume occuppied by the region.
 * @param material a pointer to a material
 */
void Region::setMaterial(Material* material) {
    _material = material;
}


/**
 * @brief Sets the squared geometric buckling for the geometry.
 * @details This method also sets the bucklking squared for the material
 *          filling it.
 * @param buckling_squared the squared geometric buckling
 */
void Region::setBucklingSquared(float buckling_squared) {
    _buckling_squared = buckling_squared;
    _material->setBucklingSquared(_buckling_squared);
}



/*******************************************************************************
 *****************************  Infinite Medium Region  ************************
 ******************************************************************************/

/**
 * @brief InfiniteMediumRegion constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_name the name of the region 
 */
InfiniteMediumRegion::InfiniteMediumRegion(const char* region_name):  
    Region(region_name) {
        _region_type = INFINITE_MEDIUM;
}


/**
 * @brief This method collides a neutron within the region.
 * @details This method encapsulates all of the neutron scattering physics
 *          which is further encapsulated by the material and isotope classes.
 * @param neutron the neutron of interest
 */
void InfiniteMediumRegion::collideNeutron(neutron* neutron) {

    if (_material == NULL)
        log_printf(ERROR, "Region %s must have material to"
			" collide neutron", _region_name);

    /* Collide the neutron in the Region's Material */
    _material->collideNeutron(neutron);
    //    neutron->_path_length = 1.0 / neutron->_total_xs;

    return;
}



/*******************************************************************************
 ************ Heterogeneous-Homogeneous Equivalence Theory Region  *************
 ******************************************************************************/

/**
 * @brief EquivalenceRegion constructor.
 * @details Sets defaults for the geometric parameters: fuel radius (0.45 cm),
 *          pitch (1.26 cm).
 * @param region_name the name of the region 
 */
EquivalenceRegion::EquivalenceRegion(const char* region_name): 
    Region(region_name) {
        _fuel_radius = 0.0;
        _pitch = 0.0;
        _half_width = 0.0;
        _other_region = NULL;
        _num_prob = 0;
}

/**
 * @brief Returns the fuel pin radius.
 * @return the fuel pin radius
 */
float EquivalenceRegion::getFuelPinRadius() {
    return _fuel_radius;
}


/**
 * @brief Returns the pin cell pitch.
 * @return the pin cell pitch
 */
float EquivalenceRegion::getPinCellPitch() {
    return _pitch;
}



/**
 * @brief This method returns the index for a certain lethargy (log10(eV)) into
 *        the uniform lethargy grid of the region's first flight
 *        collision probabilies.
 * @details The index is for the closest lethargy which is greater than or
 *          equal to the input lethargy.
 * @param lethargy the lethargy (log10(eV)) of interest
 * @return the index into the uniform lethargy grid
 */
int EquivalenceRegion::getEnergyGridIndex(float lethargy) {

    int index;

    if (lethargy > _end_lethargy)
        index = _num_prob - 1;
    else if (lethargy < _start_lethargy)
        index = 0;
    else
        index = int(floor((lethargy - _start_lethargy) / _delta_lethargy));

    return index;
}


/**
 * @brief Returns true if this region is the fuel, false otherwise.
 * @return true if fuel, false otherwise
 */
bool EquivalenceRegion::isFuel() {
    if (_region_type == EQUIVALENT_FUEL)
        return true;
    else
        return false;
}


/**
 * @brief Returns true if this region is the moderator, false otherwise.
 * @return true if moderator, false otherwise
 */
bool EquivalenceRegion::isModerator() {
    if (_region_type == EQUIVALENT_MODERATOR)
        return true;
    else
        return false;
}


/**
 * @brief Sets the first flight collision probabilities.
 * @details Sets the first flight collision (fuel-to-fuel and
 *        moderator-to-fuel) at each of the energies for which the
 *        isotope cross-sections are defined on a uniform lethargy grid.
 *        The first flight collison probabilities are computed by the
 *        geometry and set at the beginning of each PINSPEC simulation.
 * @param prob_ff an array of fuel-to-fuel collision probabilities
 * @param prob_mf an array of moderator-to-fuel collision probabilities
 * @param prob_energies an array of energies for which the first flight
 *        collision probabilities are defined
 * @param num_prob the number of first flight collision probabilities
 */
void EquivalenceRegion::setFirstFlightCollProb(float* prob_ff, float* prob_mf, 
					 float* prob_energies, int num_prob) {
    _prob_ff = prob_ff;
    _prob_mf = prob_mf;
    _prob_energies = prob_energies;
    _num_prob = num_prob;

    _start_lethargy = log10(_prob_energies[0]);
    _end_lethargy = log10(_prob_energies[_num_prob-1]);
    _delta_lethargy = (_end_lethargy - _start_lethargy) / _num_prob;
}


/**
 * @brief Sets other region for a homogeneous equivalence region.
 * @details If this region is an EQUIVALENT_FUEL region, the other region
 *          must be EQUIVALENT_MODERATOR region type. If this region is an
 *          EQUIVALENT_MODERATOR region, the other region must be
 *          EQUIVALENT_FUEL region type.
 * @param region the other region
 */
void EquivalenceRegion::setOtherRegion(EquivalenceRegion* region) {
  if (_region_type == EQUIVALENT_FUEL && 
        region->getRegionType() == EQUIVALENT_FUEL)
            log_printf(ERROR, "Unable to add an EQUIVALENT_FUEL region %s to "
	       "region %s which is also an EQUIVALENT_FUEL region type",
	       region->getName(), _region_name);
    if (_region_type == EQUIVALENT_MODERATOR && 
        region->getRegionType() == EQUIVALENT_MODERATOR)
            log_printf(ERROR, "Unable to add an EQUIVALENT_MODERATOR region "
	       "%s to region %s which is also an EQUIVALENT_MODERATOR region "
	       "type", region->getName(), _region_name);
    if (region->getRegionType() != EQUIVALENT_MODERATOR && 
        region->getRegionType() != EQUIVALENT_FUEL)
            log_printf(ERROR, "Unable to add region %s which is of %d region "
		"type to region %s since it is not a homogeneous equivalent "
		 "region", region->getName(), region->getRegionType(), 
		       _region_name);

    _other_region = region;
}


/**
 * @brief Sets the fuel radius if the region.
 * @param radius the fuel pin radius (cm)
 */
void EquivalenceRegion::setFuelPinRadius(float radius) {

    _fuel_radius = radius;

    if (_pitch != 0.0) {
        if (_region_type == EQUIVALENT_MODERATOR)
            _volume = _pitch * _pitch - M_PI * _fuel_radius * _fuel_radius;
        else
            _volume = M_PI * _fuel_radius * _fuel_radius;

        if (_material != NULL)
            _material->incrementVolume(_volume);
    }
}


/**
 * @brief Sets the pin cell pitch.
 * @param pitch the pin cell pitch (cm)
 */
void EquivalenceRegion::setPinCellPitch(float pitch) {

    _pitch = pitch;
    _half_width = pitch / 2.0;

    if (_fuel_radius != 0.0) {
        if (_region_type == EQUIVALENT_MODERATOR)
            _volume = _pitch * _pitch - M_PI * _fuel_radius * _fuel_radius;
        else
            _volume = M_PI * _fuel_radius * _fuel_radius;

        if (_material != NULL)
            _material->incrementVolume(_volume);
    }
}


/**
 * @brief Finds the first flight fuel-to-fuel collision probability
 *        for a neutron.
 * @details Uses linear interpolation to compute the first flight
 *          fuel-to-fuel collision probabilty at a neutron's energy
 *          using the region's first flight collision probability table.
 * @param neutron the neutron of interest
 * @return the first flight fuel-to-fuel collision probability
 */
float EquivalenceRegion::computeFuelFuelCollsionProb(neutron* neutron) {

    /* Find the index into the first flight collision probabilities array -
     * this index is for the nearest energy less than or equal to the
     * neutron's energy */
    float lethargy = log10(neutron->_energy);
    int lower_index = getEnergyGridIndex(lethargy);

    /* Use linear interpolation for the probability */
    float lower_prob_ff = _prob_ff[lower_index];
    float upper_prob_ff = _prob_ff[lower_index + 1];
    float delta_prob_ff = upper_prob_ff - lower_prob_ff;
    float slope = delta_prob_ff / _delta_lethargy;
    float lower_lethargy = _start_lethargy + _delta_lethargy * lower_index;
    float prob_ff =  lower_prob_ff + slope * (lethargy - lower_lethargy); 

    return prob_ff;
}


/**
 * @brief Finds the first flight moderator-to-fuel collision probability
 *        for a neutron.
 * @details Uses linear interpolation to compute the first flight
 *          moderator-to-fuel collision probabilty at a neutron's energy
 *          using the region's first flight collision probability table.
 * @param neutron the neutron of interest
 * @return the first flight moderator-to-fuel collision probability
 */
float EquivalenceRegion::computeModeratorFuelCollisionProb(neutron* neutron) {

    /* Find the index into the first flight collision probabilities array -
     * this index is for the nearest energy less than or equal to the
     * neutron's energy */
    float lethargy = log10(neutron->_energy);
    int lower_index = getEnergyGridIndex(lethargy);

    /* Use linear interpolation for the probability */
    float lower_prob_mf = _prob_mf[lower_index];
    float upper_prob_mf = _prob_mf[lower_index + 1];
    float delta_prob_mf = upper_prob_mf - lower_prob_mf;
    float slope = delta_prob_mf / _delta_lethargy;
    float lower_lethargy = _start_lethargy + _delta_lethargy * lower_index;
    float prob_mf =  lower_prob_mf + slope * (lethargy - lower_lethargy); 

    return prob_mf;
}


/**
 * @brief This method collides a neutron within the region.
 * @details This method encapsulates all of the neutron scattering physics
 *          which is further encapsulated by the material and isotope classes.
 * @param neutron the neutron of interest
 */
void EquivalenceRegion::collideNeutron(neutron* neutron) {

    float test = float(rand()) / (float)RAND_MAX;

    /* If the neutron is in the fuel */
    if (_region_type == EQUIVALENT_FUEL) {

        float prob_ff = computeFuelFuelCollsionProb(neutron);

        /* If the test is larger than prob_ff, move to moderator */
        if (test > prob_ff) {
	    neutron->_region = _other_region;
            _other_region->getMaterial()->collideNeutron(neutron);
	}
        /* Otherwise collide the neutron in the fuel's material */
        else
            _material->collideNeutron(neutron);

	//        neutron->_path_length = 1.0 / neutron->_total_xs;
    }

    /* If the neutron is in the moderator */
    else {

        float prob_mf = computeModeratorFuelCollisionProb(neutron);

        /* If the test is larger than prob_mf, move to fuel */
        if (test < prob_mf) {
	    neutron->_region = _other_region;
            _other_region->getMaterial()->collideNeutron(neutron);
	}
        /* Otherwise collide the neutron in the moderator's material */
        else
	    _material->collideNeutron(neutron);

	//        neutron->_path_length = 1.0 / neutron->_total_xs;
    }

    return;
}


/**
 * @brief EquivalenceFuelRegion constructor.
 * @param region_name the name of the region 
 */
EquivalenceFuelRegion::EquivalenceFuelRegion(const char* region_name):
    EquivalenceRegion(region_name) {
        _region_type = EQUIVALENT_FUEL;
}



/**
 * @brief EquivalenceModeratorRegion constructor.
 * @param region_name the name of the region 
 */
EquivalenceModeratorRegion::EquivalenceModeratorRegion(const char* region_name):
    EquivalenceRegion(region_name) {
        _region_type = EQUIVALENT_MODERATOR;
}


/**
 * @brief BoundedRegion constructor.
 * @param region_name the name of the region 
 */
BoundedRegion::BoundedRegion(const char* region_name): Region(region_name) { }


/**
 * @brief Adds a new halfspace of a bounding surface to a region.
 * @param halfspace an integer (-1 or +1) representing the surface's halfspace
 * @param surface a pointer to the bounding surface
 */
void BoundedRegion::addBoundingSurface(int halfspace, Surface* surface) {

    if (halfspace != -1 && halfspace != +1)
        log_printf(ERROR, "Unable to add a surface %s with halfspace %d. The "
		   "halfspace must be -1 or +1.", 
		   surface->getSurfaceName(), halfspace);

    /* Create a halfspace/surface pair and add it to the bounding surfaces
     * container for this region */
    std::pair<int, Surface*> pair=std::pair<int, Surface*>(halfspace, surface);
    _surfaces.push_back(pair);
}


/**
 * @brief Removes a halfspace of a bounding surface for a region.
 * @param halfspace an integer (-1 or +1) representing the surface's halfspace
 * @param surface a pointer to the bounding surface
 */
void BoundedRegion::removeBoundingSurface(int halfspace, Surface* surface) {

    std::pair<int, Surface*> item;
    std::vector< std::pair<int, Surface*> >::iterator test;

    /* Test whether or not the vector contains the halfspace/surface pair */
    item = std::pair<int, Surface*>(halfspace, surface);
    test = std::find(_surfaces.begin(), _surfaces.end(), item);

    /* If the vector of bounding surfaces contains the halfspace/surface
     * pair, then remove it */
    if (test != _surfaces.end())
        _surfaces.erase(test);
}


/**
 * @brief Check if this region contains some location in space.
 * @param x the x-coordinate of interest
 * @param y the y-coordinate of interest
 * @param z the z-coordinate of interest
 * @return if contained (true), otherwise (false)
 */
bool BoundedRegion::contains(float x, float y, float z) {

    int halfspace;
    Surface* surface;

    /* Loop over and query all bounding surfaces */
    std::vector< std::pair<int, Surface*> >::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        halfspace = (*iter).first;
        surface = (*iter).second;

        if (halfspace * surface->evaluate(x, y, z) < 1E-6)
            return false;
    }

    return true;
}


/**
 * @brief Check if this region contains a neutron at some location in space.
 * @param neutron the neutron of interest
 * @return if contained (true), otherwise (false)
 */
bool BoundedRegion::contains(neutron* neutron) {

    int halfspace;
    Surface* surface;

    /* Loop over and query all bounding surfaces */
    std::vector< std::pair<int, Surface*> >::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        halfspace = (*iter).first;
        surface = (*iter).second;

        if (halfspace * surface->evaluate(neutron) < -1E-6)
	    return false;
    }

    return true;
}


/**
 * @brief Checks if a neutron at some 2D location is on the region boundary.
 * @param neutron the neutron of interest
 * @return true if on the boundary, otherwise false
 */
bool BoundedRegion::onBoundary(neutron* neutron) {

    /* Loop over and query all bounding surfaces */
    std::vector< std::pair<int, Surface*> >::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
        if ((*iter).second->onSurface(neutron))
	    return true;
    }

    return false;
}



/**
 * @brief This method computes the parametrized distance along a neutron's 
 *        unit trajectory vector to the nearest bounding surface for this Region.
 * @param neutron the neutron of interest
 */
float BoundedRegion::computeParametrizedDistance(neutron* neutron) {

    float min_dist = std::numeric_limits<float>::infinity();
    float curr_dist;

    /* Loop over and query all bounding surfaces */
    std::vector< std::pair<int, Surface*> >::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        curr_dist = (*iter).second->computeParametrizedDistance(neutron);

	/* If the distance to this surface is the least found thus far,
	 * update minimum distance found */
        if (curr_dist < min_dist) {
	    min_dist = curr_dist;
	    neutron->_surface = (*iter).second;
	}
    }

    return min_dist;
}


/**
 * @brief This method collides a neutron within the region.
 * @details This method encapsulates all of the neutron scattering physics
 *          which is further encapsulated by the material and isotope classes.
 * @param neutron the neutron of interest
 */
void BoundedRegion::collideNeutron(neutron* neutron) {

    if (_material == NULL)
        log_printf(ERROR, "Region %s must have material to"
			" collide neutron", _region_name);

    float path_length = _material->sampleDistanceTraveled(neutron);
    float param_coll_dist = path_length / 
                           norm3D<float>(neutron->_u, neutron->_v, neutron->_w);
    float param_surf_dist = computeParametrizedDistance(neutron);

    /* The neutron collided within this region */
    if (param_coll_dist < param_surf_dist) {

        neutron->_region = this;
        neutron->_path_length = path_length;

        /* Update the neutron's location */
        neutron->_x += param_coll_dist * neutron->_u;
	neutron->_y += param_coll_dist * neutron->_v;
	neutron->_z += param_coll_dist * neutron->_w;
 
	/* Sample a collision type and update the neutron's energy 
	 *  and direction vector */
        _material->collideNeutron(neutron);
    }

    /* The neutron crossed a bounding surface for this region */
    else {

        /* Compute the path length to the surface */
        float delta_x = neutron->_u * param_surf_dist;
        float delta_y = neutron->_v * param_surf_dist;
	float delta_z = neutron->_w * param_surf_dist;
        path_length = norm3D<float>(delta_x, delta_y, delta_z);
        neutron->_path_length = path_length;

        /* The neutron crossed an INTERFACE type surface, so we "bump" 
	 * it across the surface with a tiny "nudge" */
        if (neutron->_surface->getBoundaryType() == INTERFACE) {
	    
            param_surf_dist += TINY_MOVE;

            /* Update the neutron's location to just beyond the surface 
	     * intersection point */
            neutron->_x += param_surf_dist * neutron->_u;
            neutron->_y += param_surf_dist * neutron->_v;
	    neutron->_z += param_surf_dist * neutron->_w;
	}

        /* The neutron crossed a REFLECTIVE boundary, so we move it the 
	 * boundary and reflect its direction of travel */
        else if (neutron->_surface->getBoundaryType() == REFLECTIVE) {

            /* Update the neutron's location to the intersection point
	     * on the surface */
            neutron->_x += param_surf_dist * neutron->_u;
            neutron->_y += param_surf_dist * neutron->_v;
	    neutron->_z += param_surf_dist * neutron->_w;

	    /* Update the neutron's trajectory vectory */
	    neutron->_surface->reflectNeutron(neutron);
        }

        /* The neutron crossed a VACUUM boundary surface so we kill it */
        else {
            
	    //FIXME: Should leakage across VACUUM boundaries affect tallies?
	    /* Kill the neutron */
	    neutron->_alive = false;
        }
    }

    return;
}


/**
 * @brief BoundedFuelRegion constructor.
 * @param region_name the name of the region 
 */
BoundedFuelRegion::BoundedFuelRegion(const char* region_name):  
    BoundedRegion(region_name) {
        _region_type = BOUNDED_FUEL;
}


/**
 * @brief Generates a series of equal area circular rings within the fuel.
 * @details This method clones the fuel region into many different equal
 *          area rings each defined by their own bounding circular surfaces.
 * @param num_rings the number of ring regions to subdivide the fuel pin into
 * 
 */
void BoundedFuelRegion::ringify(int num_rings) {
    log_printf(ERROR, "Ringify is not yet implemented for BOUNDED_FUEL type "
	       "regions.");
    return;
}


/**
 * @brief BoundedModeratorRegion constructor.
 * @param region_name the name of the region 
 */
BoundedModeratorRegion::BoundedModeratorRegion(const char* region_name):  
    BoundedRegion(region_name) {
        _region_type = BOUNDED_MODERATOR;
}


/**
 * @brief Generates a series of equal area circular rings within the moderator.
 * @details This method clones the moderator region into many different equal
 *          area rings each defined by their own bounding circular surfaces.
 * @param num_rings the number of ring regions to subdivide the moderator into
 *
 */
void BoundedModeratorRegion::ringify(int num_rings) {
    log_printf(ERROR, "Ringify is not yet implemented for BOUNDED_MODERATOR "
	       "type regions.");
    return;
}


/**
 * @brief BoundedGeneralRegion constructor.
 * @param region_name the name of the region 
 */
BoundedGeneralRegion::BoundedGeneralRegion(const char* region_name):
    BoundedRegion(region_name) {
        _region_type = BOUNDED_GENERAL;
}
