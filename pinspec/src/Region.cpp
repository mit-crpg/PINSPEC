#include "Region.h"

/**
 * @brief Region constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_name the (optional) name of the region 
 */
Region::Region(const char* region_name) {
    _region_name = (char*)region_name;
    _material = NULL;
     _volume = 1.0;
    _buckling_squared = 0.0;
}


/**
 * @brief Return the name of the region.
 * @return a character array representing this region's name
 */
char* Region::getRegionName() {
    return _region_name;
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
    _material->incrementVolume(_volume);
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
    neutron->_material = _material;
    _material->collideNeutron(neutron);

    return;
}


/**
 * @brief EquivalenceRegion constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_type the type of equivalence region (ie, EQUIVALENT_FUEL)
 * @param region_name the name of the region 
 */
EquivalenceRegion::EquivalenceRegion(regionType region_type, 
				     const char* region_name): 
    Region(region_name) {
    _region_type = region_type;
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
    }

    if (_material != NULL)
        _material->incrementVolume(_volume);
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
 * @brief This method collides a neutron within the region.
 * @details This method encapsulates all of the neutron scattering physics
 *          which is further encapsulated by the material and isotope classes.
 * @param neutron the neutron of interest
 */
void EquivalenceRegion::collideNeutron(neutron* neutron) {

    if (_material == NULL)
        log_printf(ERROR, "Region %s must have material to"
			" collide neutron", _region_name);

    /* Collide the neutron in the Region's Material */
    neutron->_material = _material;
    _material->collideNeutron(neutron);

    return;
}



/**
 * @brief BoundedRegion constructor.
 * @details Sets defaults for the geometric parameters to 0.
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

        if (halfspace * surface->evaluate(neutron) < 0)
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
 * @brief This method collides a neutron within the region.
 * @details This method encapsulates all of the neutron scattering physics
 *          which is further encapsulated by the material and isotope classes.
 * @param neutron the neutron of interest
 */
void BoundedRegion::collideNeutron(neutron* neutron) {

    if (_material == NULL)
        log_printf(ERROR, "Region %s must have material to"
			" collide neutron", _region_name);

    /* Collide the neutron in the Region's Material */
    neutron->_material = _material;
    _material->collideNeutron(neutron);

    return;
}


/**
 * @brief BoundedFuelRegion constructor.
 * @details Sets defaults for the geometric parameters to 0.
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
    return;
}


/**
 * @brief BoundedModeratorRegion constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_name the name of the region 
 */
BoundedModeratorRegion::BoundedModeratorRegion(const char* 
					       region_name):  
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
    return;
}



/**
 * @brief BoundedGeneralRegion constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_name the name of the region 
 */
BoundedGeneralRegion::BoundedGeneralRegion(const char* region_name):
    BoundedRegion(region_name) {
    _region_type = BOUNDED_GENERAL;
}
