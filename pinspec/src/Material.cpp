#include "Material.h"


int Material::_n = 1;


/**
 * @brief Material constructor.
 * @details Sets the user-defined name along with default values for the
 *          material density (0), material number density (0), material
 *          atomic mass (1), buckling (0) and volume (0).
 */
Material::Material(char* material_name) {

    int length = strlen(material_name);
    _material_name = new char[length];
   
    for (int i=0; i <= length; i++)
      _material_name[i] = material_name[i];

    _uid = _n;
    _n++;
    _material_density = 0.0;
    _material_number_density = 0.0;
    _material_atomic_mass = 1.0;
    _buckling_squared = 0.0;
    _volume = 0.0;
}


/**
 * @brief Mterial destructor.
 * @details Material does not need to delete its isotopes since SWIG
 *          handles garbage collection.
 */
Material::~Material() { 
    delete [] _material_name;
}


/**
 * @brief Returns the material name.
 * @return the name of this material
 */
char* Material::getMaterialName() {
    return _material_name;
}


/**
 * @brief Returns the unique ID auto-generated for the material.
 * @return a unique ID for the material
 */
int Material::getUid() const {
    return _uid;
}



/**
 * @brief Returns the total number density for all isotopes within the material.
 * @return the total number density \f$ \frac{at}{cm^3} \f$.
 */
float Material::getMaterialNumberDensity() {
    return _material_number_density * 1E24;
}


/**
 * @brief This method takes in a character array specifier for an isotope's
 *        name and returns a pointer to the Isotope
 * @details An example of how this would be called in python is as follows:
 *
 * @code
 *        h1 = material.getIsotope('H-1')
 * @endcode
 *
 * @param isotope the name of the isotope
 * @return a pointer to the Isotope
 */
Isotope* Material::getIsotope(char* isotope) {
    return _isotopes.at(isotope).second;
}


/**
 * @brief This method takes in a character array specifier for an isotope's
 *        name and returns a float for the Isotope's number density in 
 *        \f$ \frac{at}{cm^3} \f$.
 * @param isotope pointer to the isotope of interest
 * @return the isotope's number density
 */
float Material::getIsotopeNumDensity(Isotope* isotope) {
    return _isotopes.at(isotope->getIsotopeName()).first * 1E24;
}


/**
 * @brief This method takes in a character array specifier for an isotope's
 *        name and returns a float for the Isotope's number density in 
 *        \f$ \frac{at}{cm^3} \f$.
 * @param isotope the name of the isotope
 * @return the isotope's number density
 */
float Material::getIsotopeNumDensity(char* isotope) {
    return _isotopes.at(isotope).first * 1E24;
}


/**
 * @brief This method checks if the material contains a given isotope.
 * @param isotope a pointer to the isotope of interest
 * @return true if the material contains the isotope; otherwise false
 */
bool Material::containsIsotope(Isotope* isotope) {
    if(_isotopes.find(isotope->getIsotopeName()) == _isotopes.end())
        return false;
    else
	return true;
}


/**
 * @brief Returns the geometric buckling squared for the geometry.
 * @return the geometric buckling squared
 */
float Material::getBucklingSquared() {
    return _buckling_squared;
}


/**
 * @brief Returns the total volume of all regions containing this material.
 * @return the total volume defined by this material
 */
float Material::getVolume() {
    return _volume;
}


/**
 * @brief Returns the total number of energies for which cross-sections are
 *        defined for one of the material's cross-section types.
 * @details This method assumes that each isotope contains cross-sections
 *          defined on a uniform lethargy grid (100,000 points by default).
 *          The method queries the first isotope in the material's collection
 *          of isotopes to determine the number of cross-sections. The 
 *          method will return the number of values for 'capture', 'elastic',
 *          'fission', 'absorption' and 'total' cross-section types.
 * @param xs_type a character array for the cross-section type of interest
 * @return then number of cross-section values
 */
int Material::getNumXSEnergies(char* xs_type) {

    Isotope* isotope;

    /* If the material does not have any isotopes, throw exception */
    if (_isotopes.size() == 0) 
        log_printf(ERROR, "Unable to return the number of xs energies "
                   "for material %s since it has no isotopes", _material_name);

    isotope = _isotopes.begin()->second.second;

    return isotope->getNumXSEnergies(xs_type);
}


/**
 * @brief Fills an array with cross-section energy values.
 * @details This method is a helper function to allow PINSPEC users to
 *          get access to the material's nuclear data in Python. A user
 *          must initialize a numpy array of the correct size (ie, 
 *          a float64 array the length of the number of cross-section
 *          values) as input to this function. This function then fills
 *          the numpy array with the energy values for the isotope's
 *          cross-section data. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          num_energies = material.getNumXSEnergies()
 *          energies = numpy.zeros(num_energies)          
 *          material.retrieveXSEnergies(energies, num_energies, 'capture')
 * @endcode
 * 
 * @param energies an array to fill with the cross-section energies
 * @param num_xs the number of cross-section values
 * @param xs_type the type of cross-section
 */
void Material::retrieveXSEnergies(float* energies, int num_xs, char* xs_type) {

    Isotope* isotope;

    /* If the material does not have any isotopes, throw exception */
    if (_isotopes.size() == 0) 
        log_printf(ERROR, "Unable to return the xs energies for "
                    " material %s since it has no isotopes", _material_name);

    isotope = _isotopes.begin()->second.second;
    isotope->retrieveXSEnergies(energies, num_xs, xs_type);
}


/**
 * @brief Fills an array with macroscopic cross-section values.
 * @details This method is a helper function to allow PINSPEC users to
 *          get access to the material's nuclear data in Python. A user
 *          must initialize a numpy array of the correct size (ie, 
 *          a float64 array the length of the number of cross-section
 *          values) as input to this function. This function then fills
 *          the numpy array with the data values for one of the isotope's
 *          cross-sections. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          num_xs = material.getNumXSEnergies()
 *          xs = numpy.zeros(num_xs)          
 *          material.retrieveXS(xs, num_xs, 'capture')
 * @endcode
 * 
 * @param xs an array to fill with the macroscopic cross-section data
 * @param num_xs the number of cross-section values
 * @param xs_type the type of cross-section
 */
void Material::retrieveXS(float* xs, int num_xs, char* xs_type) {

    float* tmp_xs = new float[num_xs];

    /* If the material does not have any isotopes, throw exception */
    if (_isotopes.size() == 0) 
        log_printf(ERROR, "Unable to return a macro %s xs for material %s "
                        " since it has no isotopes", xs_type, _material_name);

    /* Initialize the macro xs to zero */
    for (int i=0; i < num_xs; i++)
        xs[i] = 0.0;

	/* Increment the cross-section type for each isotope */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter) {
	    /* Load the xs for this isotope into a temporary array */
	    iter->second.second->retrieveXS(tmp_xs, num_xs, xs_type);

            /* Add this into the macro xs for this material */
            for (int i=0; i < num_xs; i++)
                xs[i] += tmp_xs[i] * iter->second.first;
    }
}


/**
 * @brief Returns the total macroscopic cross-section for the material
 *        at some energy (eV).
 * @param energy energy of interest (eV)
 * @return the total macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getTotalMacroXS(float energy) {

    float sigma_t = 0;

    /* Increment sigma_t for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
         sigma_t += iter->second.second->getTotalXS(energy) 
	            * iter->second.first;

    return sigma_t;
}


/**
 * @brief Returns the total macroscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid 
 * @return the total macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getTotalMacroXS(int energy_index) {

    float sigma_t = 0;

    /* Increment sigma_t for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_t += iter->second.second->getTotalXS(energy_index)
	           * iter->second.first;

    return sigma_t;
}


/**
 * @brief Returns the total microscopic cross-section for the material
 *        at some energy (eV).
 * @param energy energy of interest (eV)
 * @return the total microscopic cross-section
 */
float Material::getTotalMicroXS(float energy) {

    float sigma_t = getTotalMacroXS(energy) / _material_number_density;

    return sigma_t;
}


/**
 * @brief Returns the total microscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index index into the uniform lethargy grid
 * @return the total microscopic cross-section
 */
float Material::getTotalMicroXS(int energy_index) {

    float sigma_t = 0;

    /* Increment sigma_t for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_t += iter->second.second->getTotalXS(energy_index);

    return sigma_t;
}


/**
 * @brief Returns the total macroscopic elastic scattering cross-section 
 *        for the material at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total macroscopic elastic scattering cross-section 
 *         \f$ (cm^{-1}) \f$
 * @param energy energy of interest (eV)
 * @return the total elastic macroscopic scattering cross-section (cm^-1)
 */
float Material::getElasticMacroXS(float energy) {

    float sigma_s = 0;

    /* Increment sigma_s for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_s += iter->second.second->getElasticXS(energy)
		   * iter->second.first;

    return sigma_s;
}


/**
 * @brief Returns the elastic macroscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid
 * @return the elastic macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getElasticMacroXS(int energy_index) {

    float sigma_s = 0;

    /* Increment sigma_s for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_s += iter->second.second->getElasticXS(energy_index)
		   * iter->second.first;

    return sigma_s;
}


/**
 * @brief Returns the total macroscopic elastic scattering cross-section 
 *        for the material at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total elastic microscopic scattering cross-section
 */
float Material::getElasticMicroXS(float energy) {

    float sigma_s = 0;

    /* Increment sigma_s for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_s += iter->second.second->getElasticXS(energy);

    return sigma_s;
}


/**
 * @brief Returns the elastic microscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid
 * @return the elastic microscopic cross-section
 */
float Material::getElasticMicroXS(int energy_index) {

    float sigma_s = 0;

    /* Increment sigma_s for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_s += iter->second.second->getElasticXS(energy_index);

    return sigma_s;
}


/**
 * @brief Returns the total macroscopic absorption cross-section for the 
 *        material at some energy.
 * @param energy the energy of interest (eV)
 * @return the total macroscopic absorption cross-section \f$ (cm^{-1}) \f$
 * @param energy energy of interest (eV)
 * @return the total macroscopic absorption cross-section (cm^-1)
 */
float Material::getAbsorptionMacroXS(float energy) {

    float sigma_a = 0;

    /* Increment sigma_a for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_a += iter->second.second->getAbsorptionXS(energy) *
		   iter->second.first;

    return sigma_a;
}


/**
 * @brief Returns the absorption macroscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid
 * @return the absorption macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getAbsorptionMacroXS(int energy_index) {

    float sigma_a = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_a += iter->second.second->getAbsorptionXS(energy_index)
		   * iter->second.first;

    return sigma_a;
}


/**
 * @brief Returns the total microscopic absorption cross-section for the 
          material at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total microscopic absorption cross-section
 */
float Material::getAbsorptionMicroXS(float energy) {
     float sigma_a = 0;

     /* Increment sigma_a for each isotope */
     std::map<char*, std::pair<float, Isotope*> >::iterator iter;
     for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
  	sigma_a += iter->second.second->getAbsorptionXS(energy);

     return sigma_a;
}


/**
 * @brief Returns the absorption microscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the absorption microscopic cross-section
 */
float Material::getAbsorptionMicroXS(int energy_index) {

    float sigma_a = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
 	sigma_a += iter->second.second->getAbsorptionXS(energy_index);

    return sigma_a;
}


/**
 * @brief Returns the total macroscopic capture cross-section within this 
 *        material at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total macroscopic capture cross-section \f$ (cm^{-1}) \f$
 * @param energy energy of interest (eV)
 * @return the total macroscopic capture cross-section (cm^-1)
 */
float Material::getCaptureMacroXS(float energy) {

    float sigma_c = 0;

    /* Increment sigma_a for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_c += iter->second.second->getCaptureXS(energy) *
										    iter->second.first;

    return sigma_c;
}


/**
 * @brief Returns the capture macroscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid.
 * @return the capture macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getCaptureMacroXS(int energy_index) {

    float sigma_c = 0;

    /* Increment sigma_t for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
  	sigma_c += iter->second.second->getCaptureXS(energy_index)
		   * iter->second.first;

    return sigma_c;
}


/**
 * @brief Returns the total microscopic capture cross-section for the material 
 *        at some energy (eV)
 * @param energy the energy of interest (eV)
 * @return the total microscopic capture cross-section
 */
float Material::getCaptureMicroXS(float energy) {

    float sigma_a = 0;

    /* Increment sigma_a for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
 	sigma_a += iter->second.second->getCaptureXS(energy);

    return sigma_a;
}


/**
 * @brief Returns the capture microscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the inex into the uniform lethargy grid
 * @return the capture microscopic cross-section
 */
float Material::getCaptureMicroXS(int energy_index) {

    float sigma_c = 0;

    /* Increment sigma_t for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_c += iter->second.second->getCaptureXS(energy_index);

    return sigma_c;
}


/**
 * @brief Returns the total macroscopic fission cross-section for the material
 *        at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total macroscopic fission cross-section \f$ (cm^{-1}) \f$
 * @param energy energy of interest (eV)
 * @return the total macroscopic fission cross-section (cm^-1)
 */
float Material::getFissionMacroXS(float energy) {

    float sigma_f = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
   	sigma_f += iter->second.second->getFissionXS(energy) *
		   iter->second.first;

    return sigma_f;
}


/**
 * @brief Returns the fission macroscopic cross-section within the Material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index index into the uniform lethargy grid
 * @return the fission macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getFissionMacroXS(int energy_index) {

    float sigma_f = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_f += iter->second.second->getFissionXS(energy_index)
		   * iter->second.first;

    return sigma_f;
}


/**
 * @brief Returns the total microscopic fission cross-section for the Material
 *        at some energy (eV).
 * @param energy the energy of interest (eV)
 * @return the total microscopic fission cross-section
 */
float Material::getFissionMicroXS(float energy) {

    float sigma_f = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_f += iter->second.second->getFissionXS(energy);

    return sigma_f;
}


/**
 * @brief Returns the fission microscopic cross-section for the material
 *        at some index into the uniform lethargy grid of cross-section data.
 * @param energy_index the index into the uniform lethargy grid
 * @return the fission microscopic cross-section
 */
float Material::getFissionMicroXS(int energy_index) {

    float sigma_f = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
        sigma_f += iter->second.second->getFissionXS(energy_index);

    return sigma_f;
}


/**
 * @brief Returns the total microscopic transport cross-section for the 
 *        material at some energy (eV)
 * @param energy the energy of interest (eV)
 * @return the total microscopic transport cross-section
 */
float Material::getTransportMicroXS(float energy) {

    float sigma_tr = 0;

    /* Increment sigma_a for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_tr += iter->second.second->getTransportXS(energy);

    return sigma_tr;
}


/**
 * @brief Returns the transport macroscopic cross-section for the material
 *        at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the transport macroscopic cross-section \f$ (cm^{-1}) \f$
 */
float Material::getTransportMacroXS(int energy_index) {

    float sigma_tr = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_tr += iter->second.second->getTransportXS(energy_index)
		    * iter->second.first;

    return sigma_tr;
}


/**
 * @brief Returns the total macroscopic transport cross-section for the material
 *        at some energy (eV)
 * @param energy the energy of interest (eV)
 * @return the total macroscopic transport cross-section \f$ (cm^{-1}) \f$
 */
float Material::getTransportMacroXS(float energy) {
    
    float sigma_tr = 0;

    /* Increment sigma_a for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_tr += iter->second.second->getTransportXS(energy) *
		    iter->second.first;

    return sigma_tr;
}


/**
 * @brief Returns the transport microscopic cross-section for the material
 *        at some index into the uniform lethargy grid.
 * @param energy_index the index into the uniform lethargy grid
 * @return the transport microscopic cross-section
 */
float Material::getTransportMicroXS(int energy_index) {

    float sigma_tr = 0;

    /* Increment sigma_f for each isotope */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
	sigma_tr += iter->second.second->getTransportXS(energy_index);

    return sigma_tr;
}


/**
 * @brief Sets the name for the material.
 * @param name a character array for the material's name
 */
void Material::setMaterialName(char* name) {
    _material_name = name;
}


/**
 * @brief Sets this material's density.
 * @param density the density of the material
 * @param unit the density units ('g/cc' or 'at/cc' or 'at/barncm')
 */
void Material::setDensity(float density, char* unit) {
    if (strcmp(unit, "g/cc") == 0) {
        _material_density = density;
	_density_unit = GRAM_CM3;
    }
    else if (strcmp(unit, "at/cc") == 0) {
        _material_number_density = density / 1e24;
	_density_unit = NUM_CM3;
    }
    else if (strcmp(unit, "at/barncm") == 0) {
        _material_number_density = density;
	_density_unit = NUM_BARNCM;
    }  
    else {
        log_printf(ERROR, "Cannot set Material %s number density in"
		   " units of %s since PINSPEC only support units in"
		   " g/cc, at/cc, and at/barncm", _material_name, 
		   unit);
    }    
}


/**
 * @brief Sets the material's number density.
 * @param number_density the number density of material (at/barncm)
 * @param unit the density units ('g/cc' or 'at/cc' or 'at/barncm')
 */
void Material::setNumberDensity(float density, const char* unit) 
{
    /* "g/cc" is not really a number density unit, but to be general 
     * (and to potentially save user some trouble in switching between
     * units), it is supported here anyway */
    if (strcmp(unit, "g/cc") == 0) {
        _material_density = density;
	_density_unit = GRAM_CM3;
    }
    else if (strcmp(unit, "at/cc") == 0) {
        _material_number_density = density / 1e24;
	_density_unit = NUM_CM3;
    }
    else if (strcmp(unit, "at/barncm") == 0) {
        _material_number_density = density;
	_density_unit = NUM_BARNCM;
    }  
    else {
        log_printf(ERROR, "Cannot set Material %s number density in"
		   " units %s since PINSPEC only support units in"
		   " g/cc, at/cc, and at/barncm", _material_name, unit);
    }    

}


/**
 * @brief Sets the material's total atomic mass.
 * @details If the material is a molecule this is equivalent to setting
 *          the molecular mass. For example, for UO2, the material atomic
 *           mass would be \f$ 235 + 16 + 16 = 267 \f$.
 * @param atomic_mass the material's total atomic mass
 */
void Material::setAtomicMass(float atomic_mass) {
	_material_atomic_mass = atomic_mass;
}


/**
 * @brief Sets the squared geomtric buckling for the geometry in which this
 *        material resides
 * @details This is used such that leakage can be sampled within isotopes 
 *          in the same way as all other reaction rates.
 * @param buckling_squared the squared geometric buckling for the geometry
 */
void Material::setBucklingSquared(float buckling_squared) {
    _buckling_squared = buckling_squared;
}


/**
 * @brief Increments the volume occupied by this material in the geometry.
 * @details This is used when a material is added to more than one region, as
 *          is the case for a heterogeneous pin cell with multiple radial 
 *          regions of the same material.
 * @param volume the volume to add to the total material volume
 */
void Material::incrementVolume(float volume) {
    _volume += volume;
}


/**
 * @brief Adds a new isotope to this material with a given atomic ratio.
 * @details The atomic ratio is the number of atoms of this isotope per
 *          equivalent molecule of the material. For example, for a material
 *          of UO2 one would add uranium and oxygen isotopes as follows:
 *
 * @code
 *          uo2.addIsotope(u235, 1.)
 *          uo2.addIsotope(o16, 2.)
 * @endcode
 *          
 * @param isotope a pointer to the isotope
 * @param atomic_ratio the atomic ratio of the isotope within the material
 */
void Material::addIsotope(Isotope* isotope, float atomic_ratio) {

    float isotope_number_density, N_av;
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    std::map<Isotope*, float> ::iterator iter_AO;

    /* Remove prior version of this isotope if it is already in the material */
    iter = _isotopes.find(isotope->getIsotopeName());
    if (iter != _isotopes.end())
        _isotopes.erase(iter);

    iter_AO = _isotopes_AO.find(isotope);
    if (iter_AO != _isotopes_AO.end())
        _isotopes_AO.erase(iter_AO);

    /* Add isotope and isotope AO to _isotopes_AO map */
    std::pair<Isotope*, float> new_isotope_AO =
   	std::pair<Isotope*, float> (isotope, atomic_ratio);
    _isotopes_AO.insert(new_isotope_AO);

    /* Compute the total atomic ratio */
    float total_AO = 0.0;
    for (iter_AO =_isotopes_AO.begin(); iter_AO != _isotopes_AO.end(); 
	 ++iter_AO){
    	total_AO += iter_AO->second;
    }

    /* Sum the partial contributions to the material atomic mass */
    _material_atomic_mass = 0.0;
    for (iter_AO =_isotopes_AO.begin(); iter_AO != _isotopes_AO.end(); 
	 ++iter_AO){
    	_material_atomic_mass += iter_AO->second * iter_AO->first->getA();
    }

    N_av = 6.023E-1;
    if (_density_unit == GRAM_CM3)
    {
	/* Checks to make sure material density is set already */
	if (_material_density <= 0)
	{
	    log_printf(ERROR, "Unable to add Isotope %s because the number"
		       " density for Material %s <= 0. Possible reasons: "
		       " it may not be set, or you are setting it to a"
		       " negative value", 
		       isotope->getIsotopeName(), _material_name);
	}

	/* Calculates the material's number density */
	_material_number_density = _material_density * N_av / 
	    _material_atomic_mass;
    }
    else if ((_density_unit == NUM_CM3) || (_density_unit == NUM_BARNCM))
    {
	if (_material_number_density <= 0)
	{
	    log_printf(ERROR, "Unable to add Isotope %s because the number"
		       " density for Material %s <= 0. Possible reasons: "
		       " it may not be set, or you are setting it to a"
		       " negative value", 
		       isotope->getIsotopeName(), _material_name);
	}

	_material_density = _material_number_density * _material_atomic_mass 
	    / N_av;
    }    
    else
	log_printf(ERROR, "Unable to support density unit that is not g/cc or"
		   " at/cc or at/barncm");


    /* Calculates the isotope's number density */
    isotope_number_density = atomic_ratio / total_AO * _material_number_density;

    /* Creates a pair between the number density and isotope pointer */
    std::pair<float, Isotope*> new_pair = std::pair<float, Isotope*>
	(isotope_number_density, isotope);

    std::pair<char*, std::pair<float, Isotope*> > new_isotope =
	std::pair<char*, std::pair<float, Isotope*> >
	(isotope->getIsotopeName(), new_pair);

    /* Inserts the isotope and increments the total number density */
    _isotopes.insert(new_isotope);

    /* Loop over all isotopes: update all the number densities */
    for (iter =_isotopes.begin(); iter != _isotopes.end(); ++iter){
    	/* Update isotope's number density */
    	iter->second.first = _isotopes_AO.at(iter->second.second) 
	  * _material_number_density;
        log_printf(INFO, "Isotope %s has number density %1.3E in material %s", 
                        iter->first, iter->second.first*1E24, _material_name);
    }

    return;
}


/**
 * @brief Samples a random distance to collision within this material.
 * @details Samples a random distance to the nearest collision in this
 *          material at some neutron energy using direct sampling from:
 *          \f$ \frac{-log(\xi)}{\Sigma_t(E)} \f$
 * @param neutron the neutron of interest
 */
float Material::sampleDistanceTraveled(neutron* neutron) {
    float sigma_t = getTotalMacroXS(neutron->_energy);
    return -log(float(rand()) / RAND_MAX) / sigma_t;
}


/**
 * @brief Samples an isotope for a collision.
 * @details The probability for collision with an isotope isbased on the
 *          ratios of each isotope's total cross-section to the total 
 *          cross-section of all isotope's in this Material.
 * @return a pointer to the sampled isotope 
 */
void Material::sampleIsotope(neutron* neutron) {

    float energy = neutron->_energy;
    float sigma_t = getTotalMacroXS(energy);
    //    neutron->_total_xs = sigma_t;
    neutron->_path_length = 1.0 / sigma_t;

    float sigma_t_ratio = 0.0;
    float new_sigma_t_ratio = 0.0;
    float test = float(rand()) / RAND_MAX;

    /* Loop over all isotopes */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    Isotope* isotope = NULL;

    for (iter =_isotopes.begin(); iter !=_isotopes.end(); ++iter){

        new_sigma_t_ratio += (iter->second.second->getTotalXS(energy) *
				  iter->second.first) / sigma_t;

        if (test >= sigma_t_ratio && ((test <= new_sigma_t_ratio) ||
				  fabs(test - new_sigma_t_ratio) < 1E-4)) {
            isotope = iter->second.second;
            break;
        }

        sigma_t_ratio = new_sigma_t_ratio;
    }

    if (isotope == NULL) {
        log_printf(ERROR, "Unable to find isotope type in material %s"
    		   " sampleIsotope method, energy = %1.20f, test = %1.20f,"
    		   " new_sigma_t_ratio = %1.20f", 
    		   _material_name, energy, test, new_sigma_t_ratio);
    }
    else
        neutron->_isotope = isotope;

    return;
}


/**
 * @brief Returns the material's density.
 * @return the material's density
 */
float Material::getDensity(){
    return _material_density;
}


/**
 * @brief Samples an isotope and a reaction rate for a neutron.
 * @details For a given energy, this method calls sampleIsotope() to sample 
 *        an isotope, then samples a reaction type in that isotope by using
 *        Isotope::collideNeutron() method. After this method returns, the
 *        neutron's outgoing collision energy has been updated and the neutron
 *        has been killed if it was absorbed.
 * @param neutron the neutron to collide within the material
 */
void Material::collideNeutron(neutron* neutron) {

    sampleIsotope(neutron);
    neutron->_material = this;

    neutron->_isotope->collideNeutron(neutron);

    log_printf(DEBUG, "Material %s has collided in isotope %s", 
		_material_name, neutron->_isotope->getIsotopeName());
    return;
}


/**
 * @brief This method clones a given material object by executing a deep
 *        copy of all of the material's class attributes and giving them to a 
 *        new material class object
 * @return a pointer to the new cloned material class object
 */
Material* Material::clone() {

    /* Allocate memory for the clone */
    Material* new_clone = new Material(_material_name);

    /* Set the clones atomic mass and density */
    if (_density_unit == GRAM_CM3){
        new_clone->setDensity(_material_number_density, (char*)"at/cc");
	new_clone->setDensity(_material_density, (char*)"g/cc");
    }
    else if ((_density_unit == NUM_CM3) || (_density_unit == NUM_BARNCM))
    {
        new_clone->setDensity(_material_density, (char*)"g/cc");
        new_clone->setDensity(_material_number_density, (char*)"at/cc");
    }

    /* Loops over all isotopes and add them to the clone */
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;

    for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter) {
        new_clone->addIsotope(iter->second.second, 
			      _isotopes_AO.at(iter->second.second));
    }

    new_clone->setAtomicMass(_material_atomic_mass);

    /* Return a pointer to the cloned Isotope class */
    return new_clone;
}
