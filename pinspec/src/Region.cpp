#include "Region.h"

/**
 * @brief Region constructor.
 * @details Sets defaults for the geometric parameters to 0.
 * @param region_name the name of the region 
 * @param type the type region (INFINITE, FUEL, etc)
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
 * @brief Region destructor.
 * @details The destructor does not delete anything since SWIG deals with
 *          garbage collection.
 */
Region::~Region() { }


/**
 * @brief Return the name of the region.
 * @return a character array representing this region's name
 */
char* Region::getRegionName() {
    return _region_name;
}


/**
 * @brief returns the volume of this Region \f$ (cm^3) \f$.
 * @return the region's volume
 */
float Region::getVolume() {
    return _volume;
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
 * @return the region type (INFINITE, FUEL, etc.)
 */
regionType Region::getRegionType() {
    return _region_type;
}


/**
 * @brief Returns true if this region is the fuel, false otherwise.
 * @return true if fuel, false otherwise
 */
bool Region::isFuel() {
    if (_region_type == FUEL)
        return true;
    else
        return false;
}


/**
 * @brief Returns true if this region is the moderator, false otherwise.
 * @return true if moderator, false otherwise
 */
bool Region::isModerator() {
    if (_region_type == MODERATOR)
        return true;
    else
        return false;
}


/**
 * @brief Returns true if this region is an infnite medium, false otherwise.
 * @return true if infinite, false otherwise
 */
bool Region::isInfinite() {
    if (_region_type == INFINITE)
        return true;
    else
        return false;
}


/**
 * @brief Returns the fuel pin radius if not an INFINITE region type.
 * @return the fuel pin radius
 */
float Region::getFuelRadius() {

    if (_region_type == INFINITE)
        log_printf(ERROR, "Cannot return a fuel pin radius for region %s"
		   " which is INFINITE", _region_name);

    return _fuel_radius;
}


/**
 * @brief Returns the pin cell pitch if not an INFINITE region type
 * @return the pin cell pitch
 */
float Region::getPitch() {

    if (_region_type == INFINITE)
        log_printf(ERROR, "Cannot return a pin cell pitch for region %s"
		   " which is INFINITE", _region_name);

    return _pitch;
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
 * @brief Sets the fuel radius if the region is not of INFINITE type.
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
 * @brief Sets the pin cell pitch if the region is not of INFINITE type.
 * @param pitch the pin cell pitch
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
 * @brief This method collides a neutron within the region.
 * @details This method encapsulates all of the neutron scattering physics
 *          which is further encapsulated by the material and isotope classes.
 * @param neutron the neutron of interest
 */
void Region::collideNeutron(neutron* neutron) {

    if (_material == NULL)
        log_printf(ERROR, "Region %s must have material to"
			" collide neutron", _region_name);

    /* Collide the neutron in the Region's Material */
    neutron->_material = _material;
    _material->collideNeutron(neutron);

    return;
}


/**
 * @brief Check if this region contains a neutron at some 2D location.
 * @param neutron the neutron of interest
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
 * @brief Checks if a neutron at some 2D location is on the region boundary.
 * @param neutron the neutron of interest
 * @return true if on the boundary, otherwise false
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
