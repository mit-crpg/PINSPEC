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
	_material_atomic_mass = 0.0;
	_rescaled = false;
}


/**
 * Material destructor deletes all isotopes within it
 */
Material::~Material() {

	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	Isotope* curr;

	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter) {
		curr = iter->second.second;
		delete curr;
	}
}


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
	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
		sigma_t += iter->second.second->getTotalXS(energy)
											* iter->second.first * 1E-24;

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
										* iter->second.first * 1E-24;

	return sigma_t;
}


/**
 * Returns the total microscopic cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic cross-section (barns)
 */
float Material::getTotalMicroXS(float energy) {

	float sigma_t = 0;

	/* Increment sigma_t for each isotope */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
		sigma_t += iter->second.second->getTotalXS(energy);

	return sigma_t;
}


/**
 * Returns the total microscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the total microscopic cross-section (barns)
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
												* iter->second.first * 1E-24;

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
										* iter->second.first * 1E-24;

	return sigma_e;
}


/**
 * Returns the total macroscopic elastic scattering cross-section within
 * this Material at some energy
 * @param energy energy of interest (eV)
 * @return the total elastic microscopic scattering cross-section (barns)
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
 * Returns the elastic microscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the elastic microscopic cross-section (barns)
 */
float Material::getElasticMicroXS(int energy_index) {

	float sigma_e = 0;

	/* Increment sigma_e for each isotope */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter)
		sigma_e += iter->second.second->getElasticXS(energy_index);

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
											iter->second.first * 1E-24;

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
										* iter->second.first * 1E-24;

	return sigma_a;
}


/**
 * Returns the total microscopic absorption cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic absorption cross-section (barns)
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
 * Returns the absorption microscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the absorption microscopic cross-section (barns)
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
											iter->second.first * 1E-24;

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
										* iter->second.first * 1E-24;

	return sigma_c;
}


/**
 * Returns the total microscopic capture cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic capture cross-section (barns)
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
 * Returns the capture microscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the capture microscopic cross-section (barns)
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
											iter->second.first * 1E-24;

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
										* iter->second.first * 1E-24;

	return sigma_f;
}


/**
 * Returns the total microscopic fission cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic fission cross-section (barns)
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
 * Returns the fission microscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the fission microscopic cross-section (barns)
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
 * Returns the total microscopic transport cross-section within this Material
 * at some energy
 * @param energy energy of interest (eV)
 * @return the total microscopic transport cross-section (barns)
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
										* iter->second.first * 1E-24;

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
										iter->second.first * 1E-24;

	return sigma_tr;
}


/**
 * Returns the transport microscopic cross-section within this Material
 * at some index into a rescaled energy grid
 * @param energy energy of interest (eV)
 * @return the transport microscopic cross-section (barns)
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
 * This method returns whether or not the Material's Isotope's
 * cross-sections have been rescaled to a uniform energy grid
 * @return whether or not the cross-sections have been rescaled
 */
bool Material::isRescaled() {
	return _rescaled;
}


/**
 * This method returns the index for a certain energy (eV) into
 * the uniform energy grid if this Material's Isotope's
 * cross-sections have been rescaled
 * @param energy the energy (eV) of interest
 * @return the index into the uniform energy grid
 */
int Material::getEnergyGridIndex(float energy) {

	int index;

	if (!_rescaled)
		log_printf(ERROR, "Unable to return an index for material %s "
				"since it has not been rescaled", _material_name);

	if (_scale_type == EQUAL) {
		if (energy > _end_energy)
			index = _num_energies - 1;
		else if (energy < _start_energy)
			index = 0;
		else
			index = floor((energy - _start_energy) / _delta_energy);
	}

	else if (_scale_type == LOGARITHMIC)
		energy = log10(energy);

		if (energy > _end_energy)
			index = _num_energies - 1;
		else if (energy < _start_energy)
			index = 0;
		else
			index = floor((energy - _start_energy) / _delta_energy);

	return index;
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
	if (strcmp(unit, "g/cc") != 0)
	    log_printf(ERROR, "Cannot set Material %s number density in"
						"units %s since PINSPEc only support units in" 							"g/cc", _material_name, unit);
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
    float num_density;

    /* Checks to make sure material density is set already */
    if (_material_density <= 0)
	log_printf(ERROR, "material number density is not set yet!");

    /* Increments the material's atomic mass */
    _material_atomic_mass += atomic_ratio * isotope->getA();

    /* Rescale isotope's cross sections */
    float* grid;
    grid = logspace<float, float>(_start_energy, _end_energy, _num_energies);
    isotope->rescaleXS(grid, _num_energies);
    delete [] grid;
    //_rescaled = true;

    /* Creates a pair between the number density and isotope pointer */
    std::pair<float, Isotope*> new_pair = std::pair<float, Isotope*>
	(atomic_ratio, isotope);
    
    std::pair<char*, std::pair<float, Isotope*> > new_isotope =
	std::pair<char*, std::pair<float, Isotope*> >
	(isotope->getIsotopeType(), new_pair);
    
    /* Inserts the isotope and increments the total number density */
    _isotopes.insert(new_isotope);

    return;
}


/* Update all isotopes' number densities after a material is complete */
void Material::complete() {
    std::map<char*, std::pair<float, Isotope*> >::iterator iter;
    float N_av = 6.023E-1;
    float base = _material_density * N_av / _material_number_density;
    
    /* Loop over all isotopes */
	for (iter =_isotopes.begin(); iter !=_isotopes.end(); ++iter){
	    float atomic_ratio = iter->second.first;
	    /* Computes isotope's number density */
	    float num_density = atomic_ratio * base;
	    iter->second.first = num_density;
	}

	return;
}


void Material::rescaleCrossSections(float start_energy, float end_energy,
				    int num_energies, binSpacingTypes scale_type) {

	float* grid;

	if (scale_type == EQUAL) {
		grid = linspace<float, float>(start_energy, end_energy, num_energies);
		_start_energy = start_energy;
		_end_energy = end_energy;
		_delta_energy = (_end_energy - _start_energy) / num_energies;
	}
	else {
		grid = logspace<float, float>(start_energy, end_energy, num_energies);
		_start_energy = log10(start_energy);
		_end_energy = log10(end_energy);
		_delta_energy = (_end_energy - _start_energy) / num_energies;
	}

	_num_energies = num_energies;
	_scale_type = scale_type;

	/* Loop over all isotopes */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	Isotope* isotope;

	for (iter =_isotopes.begin(); iter !=_isotopes.end(); ++iter){
		isotope = iter->second.second;
		isotope->rescaleXS(grid, num_energies);
	}

	_rescaled = true;
	delete [] grid;
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
										iter->second.first * 1E-24) / sigma_t;

		if (test >= sigma_t_ratio && ((test <= new_sigma_t_ratio) ||
							fabs(test - new_sigma_t_ratio) < 1E-5)) {
			isotope = iter->second.second;
			break;
		}
		sigma_t_ratio = new_sigma_t_ratio;
	}

	if (isotope == NULL)
		log_printf(ERROR, "Unable to find isotope type in material %s"
			   " moveNeutron method, test = %1.20f, new_num_density_ratio "
			   "= %1.20f", _material_name, test, new_sigma_t_ratio);
	
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
 * Add a tally class object to  this material's vector of Tally
 */
void Material::addTally(Tally *tally) {
    tally->setTallyDomainType(MATERIAL);
    _tallies.push_back(tally);
    return;
}

/**
 * Clear this material's vector of Tally class object pointers
 */
void Material::clearTallies() {
	_tallies.clear();
}

/**
 * For a given energy, this method calls sampleIsotope() to sample 
 * an isotop, then sample a reaction type in that isotope by using
 * Isotope::collideNeutron(), then tally the event into the appropriate
 * tally classes for that isotope if any. 
 * @param energy the incoming neutron energy (eV)
 * @return the collision type (ELASTIC, CAPTURE, FISSION)
 */
collisionType Material::collideNeutron(float energy) {
    Isotope *isotope;
    isotope = sampleIsotope(energy);
    collisionType type = isotope->getCollisionType(energy);

    /* Obtains macroscopic cross sections for this material class  */
    float total_xs = getTotalMacroXS(energy);
    float elastic_xs = getElasticMacroXS(energy);
    float absorption_xs = getAbsorptionMacroXS(energy);
    float capture_xs = getCaptureMacroXS(energy);
    float fission_xs = getFissionMacroXS(energy);
    float transport_xs = getTransportMacroXS(energy);

    /* FIXME: replace float energy with neutron struct, and update batch_num */
    int batch_num = 1;
    float sample = energy;

    /* Tallies the event into the appropriate tally classes  */
    /* Loops over all tallies and add them to the clone */
    std::vector<Tally*>::iterator iter;
	for (iter = _tallies.begin(); iter != _tallies.end(); iter ++) {
	    Tally *tally = *iter;
	    tallyType tally_type = tally->getTallyType();
	    switch (tally_type) {
	    case FLUX:
		tally->weightedTally(sample, 1.0 / total_xs, batch_num);
	    case COLLISION_RATE:
		if (type == TOTAL)
		    tally->weightedTally(sample, 1.0, batch_num);
	    case ELASTIC_RATE:
		if (type == ELASTIC)
		    tally->weightedTally(sample, elastic_xs / total_xs, 
					 batch_num);
	    case ABSORPTION_RATE:
		if (type == ABSORPTION)
		    tally->weightedTally(sample, absorption_xs / total_xs, 
					 batch_num);
	    case CAPTURE_RATE:
		if (type == CAPTURE)
		    tally->weightedTally(sample, capture_xs / total_xs, 
					 batch_num);
	    case FISSION_RATE:
		if (type == FISSION)
		    tally->weightedTally(sample, fission_xs / total_xs, 
					 batch_num);
	    case TRANSPORT_RATE:
		if (type == TRANSPORT)
		    tally->weightedTally(sample, transport_xs / total_xs, 
					 batch_num);
	    case DIFFUSION_RATE: /* FIXME */
		if (type == DIFFUSION) 
		    tally->weightedTally(sample, 
					 1.0 / (3 * transport_xs * total_xs), 
					 batch_num); 
	    case LEAKAGE_RATE:; /* FIXME */
	    }
	}



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

	/* Loops over all tallies and add them to the clone */
	for (std::vector<Tally*>::iterator it = _tallies.begin(); 
	     it != _tallies.end(); it++) {
	    new_clone->addTally(*it);
	}

	/* Loops over all isotopes and add them to the clone */
	std::map<char*, std::pair<float, Isotope*> >::iterator iter;
	Isotope* curr;

	for (iter = _isotopes.begin(); iter != _isotopes.end(); ++iter) {
	    new_clone->addIsotope(iter->second.second, iter->second.first);
	}

	/* Set the clones isotope name, atomic number, number density */
	new_clone->setMaterialName(_material_name);
	new_clone->setDensity(_material_density, (char*)"g/cc");
	new_clone->setNumberDensity(_material_number_density);
	new_clone->setAtomicMass(_material_atomic_mass);


	/* Return a pointer to the cloned Isotope class */
	return new_clone;
}

