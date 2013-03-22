/*
 * Fissioner.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Fissioner.h"

/* Default constructor for the Fissioner class
 */
Fissioner::Fissioner() {
	_num_bins = 10000;
	_E_max = 20; 	/* Default is 20 MeV */
	buildCDF();

}


/* Fissioner destructor deallocates memory for CDF
 */
Fissioner::~Fissioner() {
	if (_num_bins != 0) {
		delete [] _cdf;
		delete [] _cdf_energies;
	}
}


int Fissioner::getNumBins() {
    return _num_bins;
}

/**
 * This method sets the number of bins that we wish to use for the CDF
 * @param num_bins the number of CDF bins
 */
void Fissioner::setNumBins(int num_bins) {
	_num_bins = num_bins;
}


/**
 * This method sets the maximum energy value for the CDF
 * @param E_max the maximum CDF energy value
 */
void Fissioner::setEMax(float E_max) {
	_E_max = E_max;
}


/**
 * This method builds the CDF of the Watt Spectrum so that the
 * Fissioner can sample from it for new neutron energies
 */
void Fissioner::buildCDF() {

	/* Allocate memory to a linearly spaced array of energy values */
	_cdf_energies = linspace<float, float>(0.0, _E_max, _num_bins);

	/* Temporary array to hold the Watt spectrum values at each energy */
	float* chi = new float[_num_bins];

	/* Loop over all CDF bins and evaluate the Watt spectrum */
	for (int i=0; i < _num_bins; i++)
		chi[i] = wattSpectrum(_cdf_energies[i]);

	_cdf = new float[_num_bins];

	/* Initialize CDF by numerical integration of Watt spectrum */
	cumulativeIntegral(_cdf_energies, chi, _cdf, _num_bins, TRAPEZOIDAL);

	/* Ensure that CDF is fully normalized */
	for (int i=0; i < _num_bins; i++)
		_cdf[i] /= float(_cdf[_num_bins-1]);

	_cdf[_num_bins-1] = 1.0;

	/* Delete spectrum values */
	delete [] chi;
}


/**
 * Returns the chi value for a given energy from the Watt spectrum
 * @param energy a fission energy value (MeV)
 * @return the value of chi at that energy
 */
float Fissioner::wattSpectrum(float energy) {
	return (0.453 * exp(-1.036 * energy) * sinh(sqrt(2.29 * energy)));
}


/**
 * Sample the Fissioner's CDF and return a neutron energy from fission in MeV
 * @return a neutron energy (MeV)
 */
float Fissioner::emitNeutronMeV() {

	if (_num_bins == 0)
		log_printf(ERROR, "Unable to sample Fissioner CDF since it "
				"has not yet been created");

	return linearInterp<float, float, float>(_cdf, _cdf_energies, _num_bins,
													float(rand()) / RAND_MAX);
}


/**
 * Sample the Fissioner's CDF and return a neutron energy from fission in eV
 * @return a neutron energy (eV)
 */
float Fissioner::emitNeutroneV() {
	return emitNeutronMeV() * 1E6;
}


void Fissioner::retrieveCDF(float* cdf, int num_bins) {
    
    for (int i=0; i < _num_bins; i++)
        cdf[i] = _cdf[i];        
}


void Fissioner::retrieveCDFEnergies(float* cdf_energies, int num_bins) {
    
    for (int i=0; i < _num_bins; i++)
       cdf_energies[i] = _cdf_energies[i]; 
}
