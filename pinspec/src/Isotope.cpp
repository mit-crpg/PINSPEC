/*
 * Isotope.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */	

#include "Isotope.h"

/**
 * Isotope constructor sets default values for some isotope properties
 */
Isotope::Isotope(char* isotope_name){

	/* Default atomic number and number densities and temperature */
    _isotope_name = isotope_name;
	parseName();

	_T = 300;
	_kB = 8.617332E-5;             /* boltzmann's constant (ev / K) */
	_fissionable = false;

	/* By default this isotope has no cross-sections */
	_num_capture_xs = 0;
	_num_elastic_xs = 0;
	_num_fission_xs = 0;
	_num_absorb_xs = 0;
	_num_total_xs = 0;
	_elastic_rescaled = false;
	_capture_rescaled = false;
	_fission_rescaled = false;
    _rescaled = false;

	/* Attempt to load xs for this isotope - if the data 
	 * exists in the cross-section library */
	loadXS();	
	
    /* Rescales the isotope's cross sections */
	_start_energy = 1E-7;
	_end_energy = 2E7;
	_num_energies = 100000;
    rescaleXS(_start_energy, _end_energy, _num_energies);
 
	/* By default the thermal scattering cdfs have not been initialized */
	_use_thermal_scattering = true;
	_num_thermal_cdfs = 0;
	_num_thermal_cdf_bins = 0;

	/* FIXME: may need to update these defaults later */
	initializeThermalScattering(1E-6, 15, 1000, 15);
}



/**
 * Isotope destructor deletes array of cross section values that
 * have been assigned to this isotope
 */
Isotope::~Isotope() {
	if (_num_elastic_xs != 0) {
		delete [] _elastic_xs;
		delete [] _elastic_xs_energies;
	}
	if (_num_absorb_xs != 0) {
		delete [] _absorb_xs;
		delete [] _absorb_xs_energies;
	}
	if (_num_capture_xs != 0) {
		delete [] _capture_xs;
		delete [] _capture_xs_energies;
	}
	if (_num_fission_xs != 0) {
		delete [] _fission_xs;
		delete [] _fission_xs_energies;
	}
	if (_num_total_xs != 0) {
		delete [] _total_xs;
		delete [] _total_xs_energies;
	}
	if (_num_thermal_cdfs != 0) {
		delete [] _thermal_dist;
		for (int i=0; i < _num_thermal_cdfs; i++)
			delete [] _thermal_cdfs[i];
		delete [] _thermal_cdfs;
		delete [] _E_to_kT;
		delete [] _Eprime_to_E;
	}
}


/**
 * Parse input name and set _A value for isotope
 */
void Isotope::parseName(){

	int A = 0;
	int i = 0;
	while (_isotope_name[i] != '\0'){

		/* if - encountered, get the A value */
		if (_isotope_name[i] == '-') {

			i++;
			A = atoi(&_isotope_name[i]);
			break;
		}

		i++;
	}

	/* Set the atomic number of isotope */
	setA(A);

	log_printf(DEBUG, "Isotope %s has atomic number %i",
									_isotope_name, _A);
}



/**
 * Returns the name of the of isotope
 * @return character array with name of isotope
 */
char* Isotope::getIsotopeName() const {
	return _isotope_name;
}


/**
 * Returns the atomic number of this isotope
 * @return the atomic number
 */
int Isotope::getA() const {
    return _A;
}


/**
 * Returns the alpha ((A-1)/(A+1))^2 values for this isotope
 * @return alpha
 */
float Isotope::getAlpha() const {
    return _alpha;
}


/**
 * Returns the number density for this isotope
 * @return the number density
 */
float Isotope::getN() const {
    return _N;
}


/**
 * Returns the relative atomic amount
 * @return the relative atomic amount 
 */
float Isotope::getAO() const {
    return _AO;
}


/**
 * Return the temperature (Kelvin) for this isotope
 * @return the temperature of this isotope
 */
float Isotope::getTemperature() const {
	return _T;
}


/**
 * Return the average value of the cosine of theta for this isotope
 * in a scattering collision
 * @return the average for mu
 */
float Isotope::getMuAverage() const {
	return _mu_avg;
}


/**
 * Return whether this isotope is fissionable (true) or not (false)
 * @return whether this isotope is fissionable
 */
bool Isotope::isFissionable() const {
	return _fissionable;
}


int Isotope::getNumXSEnergies(char* xs_type) const {

    if (!strcmp(xs_type, "elastic"))
        return _num_elastic_xs;

    else if (!strcmp(xs_type, "capture"))
        return _num_capture_xs;

    else if (!strcmp(xs_type, "fission"))
        return _num_fission_xs;

    else if (!strcmp(xs_type, "absorption"))
        return _num_absorb_xs;

    else if (!strcmp(xs_type, "total"))
        return _num_total_xs;

    else 
        return 0;
}


void Isotope::retrieveXSEnergies(float* energies, int num_xs, 
                                                char* xs_type) const {

    if (!strcmp(xs_type, "elastic")) {
        for (int i=0; i < _num_elastic_xs; i++)
            energies[i] = _elastic_xs_energies[i];
    }

    if (!strcmp(xs_type, "capture")) {
        for (int i=0; i < _num_capture_xs; i++)
            energies[i] = _capture_xs_energies[i];
    }

    if (!strcmp(xs_type, "fission")) {
        for (int i=0; i < _num_fission_xs; i++)
            energies[i] = _fission_xs_energies[i];
    }

    if (!strcmp(xs_type, "absorption")) {
        for (int i=0; i < _num_absorb_xs; i++)
            energies[i] = _absorb_xs_energies[i];
    }

    if (!strcmp(xs_type, "total")) {
        for (int i=0; i < _num_total_xs; i++)
            energies[i] = _total_xs_energies[i];
    }

}


void Isotope::retrieveXS(float* xs, int num_xs, char* xs_type) const {
    
    if (!strcmp(xs_type, "elastic")) {
        for (int i=0; i < num_xs; i++)
            xs[i] = _elastic_xs[i];
    }

    if (!strcmp(xs_type, "capture")) {
        for (int i=0; i < num_xs; i++)
            xs[i] = _capture_xs[i];
    }

    if (!strcmp(xs_type, "fission")) {
        for (int i=0; i < num_xs; i++)
            xs[i] = _fission_xs[i];
    }

    if (!strcmp(xs_type, "absorption")) {
        for (int i=0; i < num_xs; i++)
            xs[i] = _absorb_xs[i];
    }

    if (!strcmp(xs_type, "total")) {
        for (int i=0; i < num_xs; i++)
            xs[i] = _total_xs[i];
    }
}


void Isotope::setMultigroupElasticXS(double* energies, int num_energies, 
                                        double* elastic_xs, int num_xs) {

	/* Check that all energies are >0 and monotonically increasing */
	double prev_energy = energies[0];

	for (int i=0 ; i < num_xs+1; i++) {

		/* Check for monotonically increasing energy */
		if (energies[i] < prev_energy)
			log_printf(ERROR, "Unable to set multigroup elastic xs for "
					"isotope %s since all xs energies must be monotonically "
											"increasing", _isotope_name);
		/* Check that energy is non-negative */		
		if (energies[i] < 0.0)
			log_printf(ERROR, "Unable to set multigroup elastic xs for "
						"isotope %s since all xs energies must be "
										"non-negative", _isotope_name);
	}

	/* Check that all xs are non-negative */
	for (int i=0; i < num_xs; i++) {
		if (elastic_xs[i] < 0.0)
			log_printf(ERROR, "Unable to set multigroup elastic xs for "
						"isotope %s since all xs values must be "
										"non-negative", _isotope_name);
	}

    log_printf(INFO, "Setting %d-group elastic xs for isotope %s", 
                                            num_xs, _isotope_name);
	
	/* Free memory for old elastic xs from ENDF */
	delete [] _elastic_xs;
	delete [] _elastic_xs_energies;

	/* The number of elastic xs values for the multigroup case */
	/* We approximate the multigroup xs as a continuous energy xs
     * using a piecewise step function */
	_num_elastic_xs = 2 * num_xs;

	/* Allocate memory for the new multigroup xs and energies */
	_elastic_xs = new float[_num_elastic_xs];
	_elastic_xs_energies = new float[_num_elastic_xs];

	/* Lowest energy xs value */
	_elastic_xs_energies[0] = energies[0];
	_elastic_xs[0] = elastic_xs[0];

	/* Create a piecewise step function to represent the multigroup xs */
	for (int i=1; i < _num_elastic_xs-1; i+=2) {
		_elastic_xs[i] = elastic_xs[i/2];
		_elastic_xs[i+1] = elastic_xs[i/2];
		_elastic_xs_energies[i] = energies[i/2] - 1E-5;
		_elastic_xs_energies[i+1] = energies[i/2] + 1E-5;
	}

	/* Highest energy xs value */
	_elastic_xs_energies[_num_elastic_xs-1] = energies[num_xs];
	_elastic_xs[_num_elastic_xs-1] = elastic_xs[num_xs-1];

	_elastic_rescaled = false;

	generateTotalXS(pow(10., _start_energy), pow(10., _end_energy), 
													_num_energies);

	return;
}


void Isotope::setMultigroupCaptureXS(double* energies, int num_energies,
                                        double* capture_xs, int num_xs) {

	/* Check that all energies are >0 and monotonically increasing */
	double prev_energy = energies[0];

	for (int i=0 ; i < num_energies; i++) {

    	/* Check for monotonically increasing energy */
		if (energies[i] < prev_energy)
			log_printf(ERROR, "Unable to set multigroup capture xs for "
					"isotope %s since all xs energies must be monotonically "
											"increasing", _isotope_name);
		/* Check that energy is non-negative */		
		if (energies[i] < 0.0)
			log_printf(ERROR, "Unable to set multigroup capture xs for "
						"isotope %s since all xs energies must be "
										"non-negative", _isotope_name);
	}

	/* Check that all xs are non-negative */
	for (int i=0; i < num_xs; i++) {
		if (capture_xs[i] < 0.0)
			log_printf(ERROR, "Unable to set multigroup capture xs for "
						"isotope %s since all xs values must be "
										"non-negative", _isotope_name);
	}


    log_printf(INFO, "Setting %d-group capture xs for isotope %s", 
                                                num_xs, _isotope_name);
	
	/* Free memory for old capture xs from ENDF */
	delete [] _capture_xs;
	delete [] _capture_xs_energies;

	/* The number of capture xs values for the multigroup case */
	/* We approximate the multigroup xs as a continuous energy xs
     * using a piecewise step function */
	_num_capture_xs = 2 * num_xs;

	/* Allocate memory for the new multigroup xs and energies */
	_capture_xs = new float[_num_capture_xs];
	_capture_xs_energies = new float[_num_capture_xs];

	/* Lowest energy xs value */
	_capture_xs_energies[0] = energies[0];
	_capture_xs[0] = capture_xs[0];

	/* Create a piecewise step function to represent the multigroup xs */
	for (int i=1; i < _num_capture_xs-1; i+=2) {
        log_printf(NORMAL, "i = %d", i);
		_capture_xs[i] = capture_xs[i/2];
		_capture_xs[i+1] = capture_xs[i/2];
		_capture_xs_energies[i] = energies[i/2] - 1E-5;
		_capture_xs_energies[i+1] = energies[i/2] + 1E-5;
	}

	/* Highest energy xs value */
	_capture_xs_energies[_num_capture_xs-1] = energies[num_xs];
	_capture_xs[_num_capture_xs-1] = capture_xs[num_xs-1];

	_capture_rescaled = false;

	generateAbsorptionXS(pow(10., _start_energy), pow(10., _end_energy), 
													_num_energies);

	generateTotalXS(pow(10., _start_energy), pow(10., _end_energy), 
													_num_energies);


	return;
}


void Isotope::setMultigroupFissionXS(double* energies, int num_energies,
                                     double* fission_xs, int num_xs) {

	/* Check that all energies are >0 and monotonically increasing */
	double prev_energy = energies[0];

	for (int i=0 ; i < num_xs+1; i++) {

		/* Check for monotonically increasing energy */
		if (energies[i] > prev_energy)
			log_printf(ERROR, "Unable to set multigroup fission xs for "
					"isotope %s since all xs energies must be monotonically "
											"increasing", _isotope_name);
		/* Check that energy is non-negative */		
		if (energies[i] < 0.0)
			log_printf(ERROR, "Unable to set multigroup fission xs for "
						"isotope %s since all xs energies must be "
										"non-negative", _isotope_name);
	}

	/* Check that all xs are non-negative */
	for (int i=0; i < num_xs; i++) {
		if (fission_xs[i] < 0.0)
			log_printf(ERROR, "Unable to set multigroup fission xs for "
						"isotope %s since all xs values must be "
										"non-negative", _isotope_name);
	}
	
    log_printf(INFO, "Setting %d-group fission xs for isotope %s", 
                                            num_xs, _isotope_name);

	/* Free memory for old capture xs from ENDF */
	delete [] _fission_xs;
	delete [] _fission_xs_energies;

	/* The number of fission xs values for the multigroup case */
	/* We approximate the multigroup xs as a continuous energy xs
     * using a piecewise step function */
	_num_fission_xs = 2 * num_xs;

	/* Allocate memory for the new multigroup xs and energies */
	_fission_xs = new float[_num_fission_xs];
	_fission_xs_energies = new float[_num_fission_xs];

	/* Lowest energy xs value */
	_fission_xs_energies[0] = energies[0];
	_fission_xs[0] = fission_xs[0];

	/* Create a piecewise step function to represent the multigroup xs */
	for (int i=1; i < _num_fission_xs-1; i+=2) {
		_fission_xs[i] = fission_xs[i/2];
		_fission_xs[i+1] = fission_xs[i/2];
		_fission_xs_energies[i] = energies[i/2] - 1E-5;
		_fission_xs_energies[i+1] = energies[i/2] + 1E-5;
	}

	/* Highest energy xs value */
	_fission_xs_energies[_num_fission_xs-1] = energies[num_xs];
	_fission_xs[_num_fission_xs-1] = fission_xs[num_xs-1];

	_fission_rescaled = false;

	generateAbsorptionXS(pow(10., _start_energy), pow(10., _end_energy), 
													_num_energies);
	generateTotalXS(pow(10., _start_energy), pow(10., _end_energy), 
													_num_energies);

	return;
}


/**
 * Returns an elastic scattering cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the elastic scattering cross-section (barns)
 */
float Isotope::getElasticXS(float energy) const {

    if (_elastic_rescaled)
        return getElasticXS(getEnergyGridIndex(energy));
	else if (_num_elastic_xs == 0)
		return 0.0;
    else
	    return linearInterp<float, float, float>(_elastic_xs_energies,
						     _elastic_xs, _num_elastic_xs, energy);
}


/**
 * Returns an elastic cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the elastic cross-section (barns)
 */
float Isotope::getElasticXS(int energy_index) const {

	if (energy_index > _num_elastic_xs) {
	    log_printf(ERROR, "Unable to retrieve elastic xs for"
		       " isotope %s since the energy index %d is out of"
		       " bounds", _isotope_name, energy_index);
	}

	return _elastic_xs[energy_index];
}


/**
 * Returns an absorption (capture plus fission) cross-section
 * value for a certain energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the absorption cross-section (barns)
 */
float Isotope::getAbsorptionXS(float energy) const {

    if (_rescaled)
        return getAbsorptionXS(getEnergyGridIndex(energy));
	else if (_num_absorb_xs == 0)
		return 0.0;
    else
	    return linearInterp<float, float, float>(_absorb_xs_energies,
								    _absorb_xs, _num_absorb_xs, energy);
}


/**
 * Returns an absorption cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the absorption cross-section (barns)
 */
float Isotope::getAbsorptionXS(int energy_index) const {

	if (energy_index > _num_absorb_xs) {
	    log_printf(ERROR, "Unable to retrieve absorption xs for"
		       " isotope %s since the energy index %d is out of"
		       " bounds", _isotope_name, energy_index);
	}

	return _absorb_xs[energy_index];
}


/**
 * Returns a capture cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the capture cross-section (barns)
 */
float Isotope::getCaptureXS(float energy) const{

    if (_capture_rescaled)
        return getCaptureXS(getEnergyGridIndex(energy));
	else if (_num_capture_xs == 0)
		return 0.0;
    else
    	return linearInterp<float, float, float>(_capture_xs_energies,
								_capture_xs, _num_capture_xs, energy);
}


/**
 * Returns a capture cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the capture cross-section (barns)
 */
float Isotope::getCaptureXS(int energy_index) const {

	if (energy_index > _num_capture_xs)
		log_printf(ERROR, "Unable to retrieve capture xs for"
				" isotope %s since the energy index %d is out of"
				" bounds", _isotope_name, energy_index);

	return _capture_xs[energy_index];
}


/**
 * Returns a fission cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the fission cross-section (barns)
 */
float Isotope::getFissionXS(float energy) const{

    if (_fission_rescaled)
        return getFissionXS(getEnergyGridIndex(energy));
	else if (_num_fission_xs == 0)
		return 0.0;
    else
    	return linearInterp<float, float, float>(_fission_xs_energies,
							_fission_xs, _num_fission_xs, energy);
}


/**
 * Returns a fission cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the fission cross-section (barns)
 */
float Isotope::getFissionXS(int energy_index) const {

	if (energy_index > _num_fission_xs) {
	    log_printf(ERROR, "Unable to retrieve fission xs for"
		       " isotope %s since the energy index %d is out of"
		       " bounds", _isotope_name, energy_index);
	}

	return _fission_xs[energy_index];
}


/**
 * Returns a total scattering cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the total cross-section (barns)
 */
float Isotope::getTotalXS(float energy) const {

	/* If the total xs has been defined explicitly, use it to
	 * linearly interpolate to find the total cross-section */
	if (_num_total_xs != 0) {
        if (_rescaled)
            return getTotalXS(getEnergyGridIndex(energy));
        else
    		return linearInterp<float, float, float>(_total_xs_energies,
									_total_xs, _num_total_xs, energy);
    }
        
	/* Otherwise loop over all xs which have been defined and
	 * add them to a total xs */
	else {

		float total_xs = 0;
		total_xs += getAbsorptionXS(energy);
		total_xs += getElasticXS(energy);
		return total_xs;
	}
}


/**
 * Returns a total scattering cross-section value for a certain index
 * into the energy array
 * @param energy the energy (eV) of interest
 * @return the total cross-section (barns)
 */
float Isotope::getTotalXS(int energy_index) const {

	/* If the total xs has been defined explicitly, use it */
	if (_num_total_xs != 0) {

	    if (energy_index > _num_total_xs) {
		log_printf(ERROR, "Unable to retrieve total xs for"
			   " isotope %s since the energy index %d is out of"
			   " bounds", _isotope_name, energy_index);
	    }

		return _total_xs[energy_index];
	}

	/* Otherwise add absorption and elastic scatter xs to get a total xs */
	else {
		float total_xs = 0;
		total_xs += getAbsorptionXS(energy_index);
		total_xs += getElasticXS(energy_index);
		return total_xs;
	}
}


/**
 * Returns the transport cross-section for this isotope at a particular energy
 * @param energy the energy (eV) of interest
 * @return the transport cross-section (barns)
 */
float Isotope::getTransportXS(float energy) const {
	return (getTotalXS(energy) - _mu_avg * getElasticXS(energy));
}


/**
 * Returns the transport cross-section for this isotope at a particular index
 * into its energy array
 * @param energy_index the index into the energy array
 * @return the transport cross-section (barns)
 */
float Isotope::getTransportXS(int energy_index) const {
	return (getTotalXS(energy_index) - _mu_avg * getElasticXS(energy_index));
}


/**
 * This method returns true if the thermal scattering distributions
 * for this isotope are to be used in the collideNeutron method
 * @return boolean if the thermal scattering distributions exist
 */
bool Isotope::usesThermalScattering() {
	return _use_thermal_scattering;
}


/**
 * This method returns whether or not the Isotope's
 * cross-sections have been rescaled to a uniform energy grid
 * @return whether or not the cross-sections have been rescaled
 */
bool Isotope::isRescaled() const {
	return _rescaled;
}


/**
 * Set the atomic number and update alpha, eta and rho
 * @param A atomic number
 */
void Isotope::setA(int A) {
    _A = A;
	_alpha = float(_A-1)/float(_A+1) * float(_A-1)/float(_A+1);
	_eta = (float(_A)+1.0) / (2.0 * sqrt(float(_A)));
	_rho = (float(_A)-1.0) / (2.0 * sqrt(float(_A)));
	_mu_avg = 2.0 / (3.0 * _A);
}


/**
 * Set the number density (at/cm^3)
 * @param N number density (at/cm^3)
 */
void Isotope::setN(float N) {
    _N = N;
}


/**
 * Set the relative atomic amount
 * @param AO relative atomic amount
 */
void Isotope::setAO(float AO) {
    _AO = AO;
}


/* Set the temperature (Kelvin)
 * @param T the temperature (Kelvin)
 */
void Isotope::setTemperature(float T) {
	_T = T;
}


void Isotope::neglectThermalScattering() {
	_use_thermal_scattering = false;
}


void Isotope::useThermalScattering() {
	_use_thermal_scattering = false;
}

/**
 * Make this a fissionable isotope
 */
void Isotope::makeFissionable() {
	_fissionable = true;
}


/**
 * Load the cross-sections from ASCII files into this isotope
 */
void Isotope::loadXS() {

    std::string directory = getXSLibDirectory();
	std::string filename;
	struct stat buffer;
	float* energies;
	float* xs_values;

	/* Set this isotope's appropriate cross-section using the data
	 * structures */


	/******************************** ELASTIC ********************************/

	/* Check whether an elastic cross-section file exists for isotope */

	filename = directory + _isotope_name + "-elastic.txt";

	if (stat(filename.c_str(), &buffer))
		log_printf(ERROR, "Unable to load elastic xs for isotope %s"
						" since no data was found in the cross-section"
						" file %s for this isotope", 
                        filename.c_str(), _isotope_name);

	log_printf(INFO, "Loading %s-elastic.txt for isotope %s", 
						_isotope_name, _isotope_name);

	/* Find the number of cross-section values in the file */
	_num_elastic_xs = getNumCrossSectionDataPoints(filename.c_str());

	/* Initialize data structures to store cross-section values */
	energies = new float[_num_elastic_xs];
	xs_values = new float[_num_elastic_xs];

	/* Parse the file into the data structures */
	parseCrossSections(filename.c_str(), energies, xs_values);

	setElasticXS(xs_values, energies, _num_elastic_xs);


	/******************************** CAPTURE ********************************/
	/* Check whether a capture cross-section file exists for isotope */

	filename = directory + _isotope_name + "-capture.txt";

	if (stat(filename.c_str(), &buffer))
		log_printf(ERROR, "Unable to load capture xs for isotope %s"
						" since no data was found in the cross-section"
						" file %s for this isotope", 
                        filename.c_str(), _isotope_name); 

	log_printf(INFO, "Loading %s-capture.txt for isotope %s", 
						_isotope_name, _isotope_name);

	/* Find the number of cross-section values in the file */
	_num_capture_xs = getNumCrossSectionDataPoints(filename.c_str());

	/* Initialize data structures to store cross-section values */
	energies = new float[_num_capture_xs];
	xs_values = new float[_num_capture_xs];

	/* Parse the file into the data structures */
	parseCrossSections(filename.c_str(), energies, xs_values);

	setCaptureXS(xs_values, energies, _num_capture_xs);

	/******************************** FISSION ********************************/
	/* Check whether a fission cross-section file exists for isotope */

	filename = directory + _isotope_name + "-fission.txt";

    /* If this isotope is fissionable and it finds it's fission xs */
	if (!stat(filename.c_str(), &buffer)) {

	log_printf(INFO, "Loading %s-fission.txt for isotope %s", 
						_isotope_name, _isotope_name);

		/* Find the number of cross-section values in the file */
		_num_fission_xs = getNumCrossSectionDataPoints(filename.c_str());

		/* Initialize data structures to store cross-section values */
		energies = new float[_num_fission_xs];
		xs_values = new float[_num_fission_xs];

		/* Parse the file into the data structures */
		parseCrossSections(filename.c_str(), energies, xs_values);

		setFissionXS(xs_values, energies, _num_fission_xs);
		makeFissionable();
	}
    
    /* If this isotope is not fissionable and it does not find a data file,
     * set the fission cross-section to zero */
    else {

		/* Initialize data structures to store cross-section values */
        _num_fission_xs = 2;
		energies = new float[_num_fission_xs];
		xs_values = new float[_num_fission_xs];

        energies[0] = 1E-7;
        energies[1] = 1E7;
        xs_values[0] = 0.0;
        xs_values[1] = 0.0;

		setFissionXS(xs_values, energies, _num_fission_xs);
    }

	return;
}


/**
 * Set the elastic cross-section for this isotope
 * @param elastic_xs a float array of microscopic elastic xs (barns)
 * @param elastic_xs_energies a float array of energies (eV)
 * @param num_elastic_xs the number of elastic xs values
 * @param type the type of angular scattering distribution
 */
void Isotope::setElasticXS(float* elastic_xs, float* elastic_xs_energies,
													int num_elastic_xs) {
    _elastic_xs = elastic_xs;
    _elastic_xs_energies = elastic_xs_energies;
    _num_elastic_xs = num_elastic_xs;
	_elastic_rescaled = false;
}


/**
 * Set the capture cross-section for this isotope
 * @param capture_xs a float array of microscopic capture xs (barns)
 * @param capture_xs_energies a float array of energies (eV)
 * @param num_capture_xs the number of capture xs values
 */
void Isotope::setCaptureXS(float* capture_xs, float* capture_xs_energies,
													int num_capture_xs) {
    _capture_xs = capture_xs;
    _capture_xs_energies = capture_xs_energies;
    _num_capture_xs = num_capture_xs;
	_capture_rescaled = false;
}


/**
 * Set the fission cross-section for this isotope
 * @param fission_xs a float array of microscopic fission xs (barns)
 * @param fission_xs_energies a float array of energies (eV)
 * @param num_fission_xs the number of fission xs values
 */
void Isotope::setFissionXS(float* fission_xs, float* fission_xs_energies,
													int num_fission_xs) {
    _fission_xs = fission_xs;
    _fission_xs_energies = fission_xs_energies;
    _num_fission_xs = num_fission_xs;
	_fission_rescaled = false;
}


void Isotope::rescaleXS(float start_energy, float end_energy,
													int num_energies) {

	float* grid;
    float* new_energies;
    float* new_xs;

	_num_energies = num_energies;
	_start_energy = start_energy;
	_end_energy = end_energy;

	grid = logspace<float, float>(start_energy, end_energy, _num_energies);

	/* Capture xs */
	if (_num_capture_xs != 0) {
		new_energies = new float[num_energies];
		memcpy(new_energies, grid, sizeof(float)*num_energies);
		new_xs = new float[num_energies];

		for (int i=0; i < num_energies; i++)
			new_xs[i] = getCaptureXS(new_energies[i]);

		_num_capture_xs = num_energies;
		delete [] _capture_xs_energies;
		delete _capture_xs;
		_capture_xs = new_xs;
		_capture_xs_energies = new_energies;
		_capture_rescaled = true;
	}

	/* Elastic xs */
	if (_num_elastic_xs != 0) {
		new_energies = new float[num_energies];
		memcpy(new_energies, grid, sizeof(float)*num_energies);
		new_xs = new float[num_energies];

		for (int i=0; i < num_energies; i++)
			new_xs[i] = getElasticXS(new_energies[i]);

		_num_elastic_xs = num_energies;
		delete [] _elastic_xs_energies;
		delete _elastic_xs;
		_elastic_xs = new_xs;
		_elastic_xs_energies = new_energies;
		_elastic_rescaled = true;
	}

	/* Fission xs */
	if (_num_fission_xs != 0) {
		new_energies = new float[num_energies];
		memcpy(new_energies, grid, sizeof(float)*num_energies);
		new_xs = new float[num_energies];


		for (int i=0; i < num_energies; i++)
			new_xs[i] = getFissionXS(new_energies[i]);

		_num_fission_xs = num_energies;
		delete [] _fission_xs_energies;
		delete _fission_xs;
		_fission_xs = new_xs;
		_fission_xs_energies = new_energies;
		_fission_rescaled = true;
	}

	_start_energy = log10(start_energy);
	_end_energy = log10(end_energy);
	_delta_energy = (_end_energy - _start_energy) / _num_energies;

	generateAbsorptionXS(start_energy, end_energy, num_energies);
	generateTotalXS(start_energy, end_energy, num_energies);

	delete [] grid;

	return;
}


void Isotope::generateAbsorptionXS(float start_energy, float end_energy, 
														int num_energies) {

	/* Generate an absorption xs on the uniform energy/lethargy grid */
	float* new_energies = logspace<float, float>(start_energy, 
												end_energy, num_energies);
	float* new_xs = new float[num_energies];

	_num_absorb_xs = num_energies;

    if (_fissionable)
	    for (int i=0; i < num_energies; i++)
		    new_xs[i] = getCaptureXS(new_energies[i]) + 
                                                getFissionXS(new_energies[i]);

    else
	    for (int i=0; i < num_energies; i++)
		    new_xs[i] = getCaptureXS(new_energies[i]);


	_absorb_xs = new_xs;
	_absorb_xs_energies = new_energies;

	return;
}


void Isotope::generateTotalXS(float start_energy, float end_energy, 
														int num_energies) {

	/* Generate a total xs on the uniform energy/lethargy grid */
	float* new_energies = logspace<float, float>(start_energy, 
										end_energy, num_energies);
	float* new_xs = new float[num_energies];

	_num_total_xs = num_energies;

	for (int i=0; i < num_energies; i++)
		new_xs[i] = getAbsorptionXS(new_energies[i]) + 
                                                getElasticXS(new_energies[i]);

	_total_xs = new_xs;
	_total_xs_energies = new_energies;

	_rescaled = true;
}



/**
 * This method clones a given Isotope class object by executing a deep
 * copy of all of the Isotope's class attributes and giving them to a new
 * Isotope class object
 * @return a pointer to the new cloned Isotope class object
 */
Isotope* Isotope::clone() {

	/* Allocate memory for the clone */
	Isotope* new_clone = new Isotope(_isotope_name);

	/* Set the clone's atomic number, number density */
	new_clone->setA(_A);
	new_clone->setN(_N);
	new_clone->setTemperature(_T);

	/* Return a pointer to the cloned Isotope class */
	return new_clone;
}


/**
 * For a given energy, this method determines a random collision type
 * based on the values of each of its cross-section types at that energy
 * @param neutron a poiner to a nuetron
 */
void Isotope::sampleCollisionType(neutron* neutron) {

	float energy = neutron->_energy;
	float test = float(rand()) / RAND_MAX;
	float collision_xs = 0.0;
	float prev_collision_xs = 0.0;
	float total_xs = getTotalXS(energy);

    /* Elastic scatter collision */
    collision_xs += getElasticXS(energy) / total_xs;
    if (test >= prev_collision_xs && test <= collision_xs)
		return;

    /* Capture collision */
    prev_collision_xs = collision_xs;
    collision_xs += getCaptureXS(energy) / total_xs;    
    if (test >= prev_collision_xs && test <= collision_xs) {
    	neutron->_alive = false;
		return;
	}

    /* Otherise, return a fission collision */
   	neutron->_alive = false;
	return;
}


/**
 * For a given neutron energy in eV in a scattering collision, this
 * function returns the outgoing energy in eV, Eprime, for the collision
 * based on the thermal scattering distributions
 * @param energy the incoming energy (eV)
 * @return the outgoing energy (eV)
 */
float Isotope::getThermalScatteringEnergy(float energy) {

	/* First check that the thermal scattering CDFs have been initialized */
    if (_num_thermal_cdfs == 0) {
	log_printf(ERROR, "Unable to sample the thermal scattering CDFs for"
		   " isotope %s because they have not yet been initialized",
		   _isotope_name);
    }

	/* Convert energies in eV to eV / kT */
	energy /= (_kB * _T);

	/* Compute possible values for E to scatter to */
	float* possible_Eprimes = new float[_num_thermal_cdf_bins];
	for (int i=0; i < _num_thermal_cdf_bins; i++)
		possible_Eprimes[i] = _Eprime_to_E[i] * energy;

	float rn = float(rand()) / RAND_MAX;
	int index;
	float Eprime;

	/* Check if energy is lower than all thermal scattering CDFs */
	if (energy < _E_to_kT[0]) {
		index = findUpperIndex(_thermal_cdfs[0],
								_num_thermal_cdf_bins-1, 0, rn);
		Eprime = possible_Eprimes[index];
	}

	/* Check if energy is above all thermal scattering CDFs */
	else if (energy > _E_to_kT[_num_thermal_cdfs-1]) {
		index = findUpperIndex(_thermal_cdfs[_num_thermal_cdfs-1],
									_num_thermal_cdf_bins-1, 0, rn);
		Eprime = possible_Eprimes[index];
	}

	/* Otherwise the energy is sandwiched within the scattering CDFs */
	else {
		int upper_index = findUpperIndex(_E_to_kT, _num_thermal_cdfs-1,
															0, energy);
		int lower_index = upper_index - 1;
		int Eprime_lower_index = findUpperIndex(_thermal_cdfs[lower_index],
										_num_thermal_cdf_bins-1, 0, rn);
		float Eprime_lower = possible_Eprimes[Eprime_lower_index];
		int Eprime_upper_index = findUpperIndex(_thermal_cdfs[upper_index],
										_num_thermal_cdf_bins-1, 0, rn);
		float Eprime_upper = possible_Eprimes[Eprime_upper_index];
		float delta_E_to_kT = _E_to_kT[upper_index] - _E_to_kT[lower_index];
		float delta_Eprime = Eprime_upper - Eprime_lower;
		float slope = delta_Eprime / delta_E_to_kT;
		Eprime = slope * (energy - _E_to_kT[lower_index]) + Eprime_lower;
	}

	/* Convert outgoing energy back into eV */
	Eprime *= (_kB * _T);

	delete [] possible_Eprimes;

	return Eprime;
}


/* This method initializes the probability distributions for thermal
 * scattering. It takes in arguments for the start_energy and end
 * energy (ratios of kT) and the number of distributions which it
 * uses to generate logarithmically spaced energies for the distributions.
 * It also takes in the number of energy bins to use for each bin.
 * @param start_energy the first distribution's energy
 * @param end_energy the final distribution's energy
 * @param num_bins the number of bins per distribution
 * @param end_distributions the number of scattering distributions
 */
void Isotope::initializeThermalScattering(float start_energy,
					float end_energy, int num_bins, 
					  int num_distributions) {

	/* Number of thermal scattering distributions */
	_num_thermal_cdfs = num_distributions;

	/* Number of bins per distribution */
	_num_thermal_cdf_bins = num_bins;

	/* Allocate memory for distributions */
	_thermal_cdfs = new float*[_num_thermal_cdfs];
	for (int i=0; i < _num_thermal_cdfs; i++)
		_thermal_cdfs[i] = new float[_num_thermal_cdf_bins];

	_thermal_dist = new float[_num_thermal_cdfs * _num_thermal_cdf_bins];
	float* cdf = new float[_num_thermal_cdf_bins];

	/* Initialize logarithmically spaced E/kT for each distribution */
	_E_to_kT = logspace<float, float>(start_energy/(_kB*_T),
					  end_energy/(_kB*_T), 
					  _num_thermal_cdfs);

	/* Find the maximum Eprime / E value that we must extend our distributions
	 * to before they all fall below some tolerance */
	bool tolerance_met = false;
	float dist_tolerance = 0.1;
	float curr_prob;
	float curr_Eprime_to_E = 1.0;
	while (!tolerance_met) {

		/* Start with the tolerance being met */
		tolerance_met = true;

		/* Loop over all CDFs and check if we are within the threshold */
		for (int i=0; i < _num_thermal_cdfs; i++) {
			curr_prob = thermalScatteringProb(curr_Eprime_to_E, i);

			/* If we are above the tolerance */
			if (curr_prob > dist_tolerance)
				tolerance_met = false;
		}

		/* Update distance along x-axis */
		if(!tolerance_met)
			curr_Eprime_to_E += 0.25;
	}

	/* Initialize x-axis of Eprime to E ratios */
	_Eprime_to_E = logspace<float, float>(1E-5, curr_Eprime_to_E,
					      _num_thermal_cdf_bins);

	/* Loop over each distribution */
	for (int i=0; i < _num_thermal_cdfs; i++) {
		for (int j=0; j < _num_thermal_cdf_bins; j++)
			_thermal_dist[i*_num_thermal_cdf_bins + j] =
			    thermalScatteringProb(_Eprime_to_E[j], i);
	}

	/* Create CDFs for each distribution */
	for (int i=0; i < _num_thermal_cdfs; i++) {
		cumulativeIntegral(_Eprime_to_E,
				   &_thermal_dist[i*_num_thermal_cdf_bins], cdf,
				   _num_thermal_cdf_bins, TRAPEZOIDAL);

		/* Transfer CDF values to our array */
		for (int j=0; j < _num_thermal_cdf_bins; j++)
			_thermal_cdfs[i][j] = cdf[j];
	}


	delete [] cdf;

	/* Normalize CDFs */
	for (int i=0; i < _num_thermal_cdfs; i++) {
		for (int j=0; j < _num_thermal_cdf_bins; j++)
			_thermal_cdfs[i][j] /= _thermal_cdfs[i][_num_thermal_cdf_bins-1];
	}

	return;
}


/**
 * This function computes the thermal scattering probability for
 * for a ratio of initial to final energies
 * @param E_prime_to_E a ratio of initial to final energies
 * @return the probability of the ratio occurring
 */
float Isotope::thermalScatteringProb(float E_prime_to_E, int dist_index) {

	double prob;

    /* Computes the final energy for each of the ratios */
    float Eprime = _E_to_kT[dist_index] * E_prime_to_E;

    /* Uses the equation from 22.211 slide 26 of the 2nd lecture
     * to compute probabilities */
    double a = sqrt(_E_to_kT[dist_index]);
	double b = sqrt(Eprime);
	double c = erf(_eta * b - _rho * a);
	double d = erf(_eta * b + _rho * a);
	double e = erf(_eta * a - _rho * b);
	double f = erf(_eta * a + _rho * b);
	double g = exp(double(_E_to_kT[dist_index]) - double(Eprime));

	/* Account for lower and upper signs in equation */
	if (Eprime > _E_to_kT[dist_index])
		prob = (c - d) + g * (e + f);
	else
		prob = (c + d) + g * (e - f);

	/* Multiply by eta / 2 */
	prob *= double(_eta*_eta) / 2.0;

	/* Normalize to the atomic mass by multiplying by 1-alpha */
	prob *= (1.0 - double(_alpha));

	return float(prob);
}


int Isotope::getNumThermalCDFs() {
    return _num_thermal_cdfs;
}


int Isotope::getNumThermalCDFBins() {
    return _num_thermal_cdf_bins;
}


void Isotope::retrieveThermalCDFs(float* cdfs, int num_values) {
  
    for (int i=0; i < _num_thermal_cdfs; i++) {
        for (int j=0; j < _num_thermal_cdf_bins; j++)
            cdfs[i*_num_thermal_cdf_bins + j] = _thermal_cdfs[i][j];
    }  
}

void Isotope::retrieveThermalDistributions(float* dist, int num_values) {

    for (int i=0; i < _num_thermal_cdfs; i++) {
        for (int j=0; j < _num_thermal_cdf_bins; j++)
            dist[i*_num_thermal_cdf_bins + j] = _thermal_dist[i*_num_thermal_cdf_bins + j];
    }  
}



void Isotope::retrieveEtokT(float* E_to_kT, int num_cdfs) {
    
    for (int i=0; i < _num_thermal_cdfs; i++)
        E_to_kT[i] = _E_to_kT[i];
}


void Isotope::retrieveEprimeToE(float* Eprime_to_E, int num_bins) {

    for (int i=0; i < _num_thermal_cdf_bins; i++)
        Eprime_to_E[i] = _Eprime_to_E[i];
}



/**
 * For a given neutron, this method samples a distance for the neutron.
 * @param neutron structure
 * @return the distance traveled
 */
float Isotope::getDistanceTraveled(neutron* neutron) {

    double sigma_a;
	double random;
	double distance;

    sigma_a = getTotalXS(neutron->_energy);
    random = (float)(rand()) / (float)(RAND_MAX);
    distance = - log(random) / sigma_a;

    return distance;
}


/**
 * For a given energy, this method calls getCollisionType() to sample
 * the collision type, and then tally the event into the appropriate
 * tally classes for that isotope if any.
 * @param neutron structure
 */
void Isotope::collideNeutron(neutron* neutron) {

	neutron->_old_energy = neutron->_energy;

    /* obtain collision type */
    sampleCollisionType(neutron);

    /* Sample outgoing energy uniformally between [alpha*E, E] */
    float alpha = getAlpha();
    double random = (float)(rand()) / (float)(RAND_MAX);

    /* Asymptotic elastic scattering above 4 eV or if no thermal scattering */
    if (neutron->_energy > 4.0 || !_use_thermal_scattering)
    	neutron->_energy *= (alpha + (1.0 - alpha) * random);
    else
    	neutron->_energy = getThermalScatteringEnergy(neutron->_energy);

	return;
}



