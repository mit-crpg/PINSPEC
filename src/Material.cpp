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
float Material::getTotalNumberDensity() {
	return _tot_num_density;
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
 * @param name the name of this Material
 */
void Material::setMaterialName(char* name) {
	_material_name = name;
}


/**
 * Adds a new isotope to this Material
 * @param isotope a pointer to a isotope class object
 * @param num_density the number density (at/cm^3) for this isotope
 */
void Material::addIsotope(Isotope* isotope, float num_density) {

	/* Creates a pair between the number density and isotope pointer */
	std::pair<float, Isotope*> new_pair = std::pair<float, Isotope*>
													(num_density, isotope);

	std::pair<char*, std::pair<float, Isotope*> > new_isotope =
							std::pair<char*, std::pair<float, Isotope*> >
									(isotope->getIsotopeType(), new_pair);

	/* Inserts the isotope and increments the total number density */
	_isotopes.insert(new_isotope);
	_tot_num_density += num_density;

	return;
}



void Material::rescaleCrossSections(float start_energy, float end_energy,
								int num_energies, energyGridType scale_type) {

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
