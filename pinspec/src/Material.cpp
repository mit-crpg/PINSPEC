/*
 * Material.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Material.h"


/**
 * Material constructor sets empty default material name
 */
Material::Material() {
	_material_name = (char*)"";
	_material_density = 0.0;
	_material_number_density = 0.0;
	_material_atomic_mass = 1.0;
}


/**
 * Material destructor does not delete anything since that is left to SWIG
 */
Material::~Material() { }


/**
 * Returns the name of this Material as specified by the user
 * @return the name of this Material
 */
char* Material::getMaterialName() {
	return _material_name;
}


/**
 * Returns the total number density for all isotopes within
 * this material
 * @return the total number density (at/cm^3)
 */
float Material::getMaterialNumberDensity() {
	return _material_number_density;
}


/**
 * This method takes in a character array specifier for an Isotope's
 * name and returns a pointer to the Isotope
 * @param isotope the name of the isotope
 * @return a pointer to the Isotope
 */
Isotope* Material::getIsotope(char* isotope) {
	return _isotopes.at(isotope).second;
}


/**
 * This method takes in a character array specifier for an Isotope's
 * name and returns a float for the Isotope's number density in at/cm^3
 * @param isotope the name of hte isotope
 * @return the isotope's number density
 */
float Material::getIsotopeNumDensity(char* isotope) {
	return _isotopes.at(isotope).first;
}


int Material::getNumXSEnergies() const {

    Isotope* isotope;

    if (_isotopes.size() == 0) 
        log_printf(ERROR, "Unable to return the number of xs energies for "
                    " material %s since it has no isotopes", _material_name);

    isotope = _isotopes.begin()->second.second;
    return isotope->getNumXSEnergies();
}


binSpacingType Material::getEnergyGridScaleType() const {
    
    if (_isotopes.size() == 0)
        log_printf(ERROR, "Unable to get the energy grid scale type for "
                        "Material %s since it does not contain any "
                        "isotopes yet", _material_name);

    Isotope* isotope = _isotopes.begin()->second.second;
    return isotope->getEnergyGridScaleType();
}


void Material::retrieveXSEnergies(float* energies, int num_xs) const {

    Isotope* isotope;

    if (_isotopes.size() == 0) 
        log_printf(ERROR, "Unable to return the xs energies for "
                    " material %s since it has no isotopes", _material_name);

    isotope = _isotopes.begin()->second.second;
    isotope->retrieveXSEnergies(energies, num_xs);
}


void Material::retrieveXS(float* xs, int num_xs, char* xs_type) {

    float* tmp_xs = new float[num_xs];

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
 * Returns the total macroscopic cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic cross-section (cm^-1)
 */
float Material::getTotalMacroXS(float energy) {

	float sigma_t = 0;

	/* Increment sigma_t for each isotope */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter) {
		sigma_t += iter->second.second->getTotalXS(energy)
											* iter->second.first;
    }

	return sigma_t;
}


/**
 * Returns the total macroscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the total macroscopic cross-section (cm^-1)
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
 * Returns the total macroscopic elastic scattering cross-section within
 * this Material at some energy
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
 * Returns the elastic macroscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the elastic macroscopic cross-section (cm^-1)
 */
float Material::getElasticMacroXS(int energy_index) {

	float sigma_e = 0;

	/* Increment sigma_e for each isotope */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
		sigma_e += iter->second.second->getElasticXS(energy_index)
										* iter->second.first;

	return sigma_e;
}


/**
 * Returns the total macroscopic absorption cross-section within this Material
 * at some energy
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
 * Returns the absorption macroscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the absorption macroscopic cross-section (cm^-1)
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
 * Returns the total macroscopic capture cross-section within this Material
 * at some energy
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
 * Returns the capture macroscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the capture macroscopic cross-section (cm^-1)
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
 * Returns the total macroscopic fission cross-section within this Material
 * at some energy
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
 * Returns the fission macroscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the fission macroscopic cross-section (cm^-1)
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
 * Returns the transport macroscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the transport macroscopic cross-section (cm^-1)
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
 * Returns the total macroscopic transport cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total macroscopic transport cross-section (cm^-1)
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
 * Sets this Material's name as defined by the user
 * @param set the name of this Material
 */
void Material::setMaterialName(char* name) {
	_material_name = name;
}


/**
 * Sets this Material's density as defined by the user
 * @param set the density of this Material
 */
void Material::setDensity(float density, char* unit) {
	_material_density = density;
	if (strcmp(unit, "g/cc") != 0) {
	    log_printf(ERROR, "Cannot set Material %s number density in"
		       "units %s since PINSPEc only support units in"
		       "g/cc", _material_name, unit);
	}    
}


/**
 * Set the number density of this material.
 */
void Material::setNumberDensity(float number_density) {
	_material_number_density = number_density;
}


/**
 * Sets this Material's atomic mass computed from user inputs
 * @param set the atomic mass of this Material
 */
void Material::setAtomicMass(float atomic_mass) {
	_material_atomic_mass = atomic_mass;
}


/**
 * Adds a new isotope to this Material
 * @param isotope a pointer to a isotope class object
 * @param atomic_ratio the atomic ratio of the isotope
 */
void Material::addIsotope(Isotope* isotope, float atomic_ratio) {

	float isotope_number_density, N_av;
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    std::map<Isotope*, float> ::iterator iter_AO;

    /* Add isotope and isotope AO to _isotopes_AO map */
    std::pair<Isotope*, float> new_isotope_AO =
   	std::pair<Isotope*, float> (isotope, atomic_ratio);
    _isotopes_AO.insert(new_isotope_AO);

    /* Checks to make sure material density is set already */
    if (_material_density <= 0)
	log_printf(ERROR, "Unable to add Isotope %s since the number density "
                       "for Material %s has not yet been set", 
                        isotope->getIsotopeName(), _material_name);

    /* Increments the material's total atomic mass and number density */
    N_av = 6.023E-1;

    /* Compute the total atomic ratio */
    float total_AO = 0.0;
    for (iter_AO =_isotopes_AO.begin(); iter_AO != _isotopes_AO.end(); ++iter_AO){
    	total_AO += iter_AO->second;
    }

    /* Sum the partial contributions to the material atomic mass */
    _material_atomic_mass = 0.0;
    for (iter_AO =_isotopes_AO.begin(); iter_AO != _isotopes_AO.end(); ++iter_AO){
    	_material_atomic_mass += iter_AO->second / total_AO * iter_AO->first->getA();
    }

    /* Calculates the material's number density */
    /* Notice I am using old_atomic_mass because I update all isotopes at
     * the end of this function. */
    _material_number_density = _material_density * N_av / _material_atomic_mass;

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
    	iter->second.first = _isotopes_AO.at(iter->second.second) / total_AO * _material_number_density;
    }

    return;
}



/**
 * Samples a isotope for a collision with a probability based on the
 * ratios of each isotope's total cross-section to the total cross-section
 * of all isotope's in this Material
 * @return a pointer to the chosen isotope
 */
Isotope* Material::sampleIsotope(float energy) {

	float sigma_t = getTotalMacroXS(energy);
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
							fabs(test - new_sigma_t_ratio) < 1E-5)) {
			isotope = iter->second.second;
			break;
		}
		sigma_t_ratio = new_sigma_t_ratio;
	}

	if (isotope == NULL) {
	    log_printf(ERROR, "Unable to find isotope type in material %s"
		       " sampleIsotope method, test = %1.20f," 
		       " new_num_density = %1.20f", 
		       _material_name, test, new_sigma_t_ratio);
	}

	return isotope;
}


/**
 * Get the material density
 * @return the material density
 */
float Material::getDensity(){
	return _material_density;
}


/**
 * For a given energy, this method calls sampleIsotope() to sample 
 * an isotop, then sample a reaction type in that isotope by using
 * Isotope::collideNeutron(), then tally the event into the appropriate
 * tally classes for that isotope if any. 
 * @param energy the incoming neutron energy (eV)
 * @return the collision type (ELASTIC, CAPTURE, FISSION)
 */
collisionType Material::collideNeutron(neutron* neut) {

    float sample = neut->_energy;

    Isotope *isotope;
    isotope = sampleIsotope(sample);
    collisionType type = isotope->collideNeutron(neut);
    log_printf(DEBUG, "Material %s has sampled collision type %d from ",
	       "isotope %s", _material_name, type, isotope->getIsotopeName());

	return type;
}


/**
 * This method clones a given Material class object by executing a deep
 * copy of all of the Material's class attributes and giving them to a new
 * Material class object
 * @return a pointer to the new cloned Material class object
 */
Material* Material::clone() {

	/* Allocate memory for the clone */
	Material* new_clone = new Material();

	/* Loops over all isotopes and add them to the clone */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;

	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter) {
	    new_clone->addIsotope((iter->second.second)->clone(), iter->second.first);
	}

	/* Set the clones isotope name, atomic number, number density */
	new_clone->setMaterialName(_material_name);

    if (_density_unit == GRAM_CM3)
	    new_clone->setDensity(_material_density, (char*)"g/cc");
    else if (_density_unit == NUM_CM3)
	    new_clone->setDensity(_material_density, (char*)"at/cc");

	new_clone->setNumberDensity(_material_number_density);
	new_clone->setAtomicMass(_material_atomic_mass);


	/* Return a pointer to the cloned Isotope class */
	return new_clone;
}


