/*
 * Fissioner.h
 *
 *  Created on: Mar 21, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef FISSIONER_H_
#define FISSIONER_H_

#include "log.h"
#include "integrate.h"
#include "interpolate.h"
#include "arraycreator.h"


/* Fissioner - This class represents the physics of neutron emissions
 * The Fissioner class generates a cumulative distribution function
 * for the chi spectrum of fission neutron energies given by the Watt
 * spectrum given in the "Fundamentals of Nuclear Reactor Physics",
 * E. E. Lewis. This class can then use the cdf to sample neutron
 * emission energies.
 */
class Fissioner {

private:
	int _num_bins;
	float* _cdf;
	float* _cdf_energies;
	float _E_max;

public:
	Fissioner();
	virtual ~Fissioner();
	void setNumBins(int num_bins);
	void setEMax(float E_max);
	void buildCDF();
	float wattSpectrum(float energy);
	float emitNeutronMeV();
	float emitNeutroneV();
};

#endif /* FISSIONER_H_ */
