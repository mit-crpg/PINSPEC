/*
 * Isotope.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef ISOTOPE_H_
#define ISOTOPE_H_

#ifdef __cplusplus
#include <limits>
#include <vector>
#include <map>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "interpolate.h"
#include "integrate.h"
#include "arraycreator.h"
#include "xsreader.h"
#include "log.h"
#include "Neutron.h"


/**
 * The isotope class represents an isotope and all of its properties
 * which are relevant to neutronics
 */
class Isotope {
	
private:
	char* _isotope_name;
	int _A;
	float _alpha;
	float _eta;
	float _rho;
	float _N;
	float _AO;
	float _T;
	float _mu_avg;
	bool _fissionable;
	bool _rescaled;

	int _num_elastic_xs;
	float* _elastic_xs;
	float* _elastic_xs_energies;
	bool _elastic_rescaled;

	int _num_capture_xs;
	float* _capture_xs;
	float* _capture_xs_energies;
	bool _capture_rescaled;

	int _num_fission_xs;
	float* _fission_xs;
	float* _fission_xs_energies;
	bool _fission_rescaled;

	int _num_absorb_xs;
	float* _absorb_xs;
	float* _absorb_xs_energies;

	int _num_total_xs;
	float* _total_xs;
	float* _total_xs_energies;

	/* Values related to rescaled cross-sections on a uniform lethargy grid */
	int _num_energies;
	float _start_energy;
	float _end_energy;
	float _delta_energy;

	bool _use_thermal_scattering;
	int _num_thermal_cdfs;
	int _num_thermal_cdf_bins;
	float* _thermal_dist;
	float** _thermal_cdfs;
	float* _E_to_kT;
	float* _Eprime_to_E;
	float _kB;

	void loadXS();
    void loadXS(char* xs_type);
	void setElasticXS(float* elastic_xs, float* elastic_xs_energies,
			  									int num_elastic_xs);
	void setCaptureXS(float* capture_xs, float* capture_xs_energies,
											     int num_capture_xs);
	void setFissionXS(float* fission_xs, float* fission_xs_energies,
											 	 int num_fission_xs);
	void rescaleXS(float start_energy, float end_energy, int num_energies);
	void generateAbsorptionXS(float start_energy, float end_energy, int num_energies);
	void generateTotalXS(float start_energy, float end_energy, int num_energies);

	void initializeThermalScattering(float start_energy, float end_energy,
					 int num_bins, int num_distributions);
	float thermalScatteringProb(float E_prime_to_E, int dist_index);

public:
	Isotope(char *_isotope_name);
	virtual ~Isotope();

	void parseName();
	void makeFissionable();

	char* getIsotopeName() const;
	int getA() const;
    float getAlpha() const;
    float getN() const;
    float getAO() const;
    float getTemperature() const;
    float getMuAverage() const;
	bool isFissionable() const;

    int getNumXSEnergies(char* xs_type) const;

    float getElasticXS(float energy) const;
    float getElasticXS(int energy_index) const;
    float getAbsorptionXS(float energy) const;
    float getAbsorptionXS(int energy_index) const;
    float getCaptureXS(float energy) const;
    float getCaptureXS(int energy_index) const;
    float getFissionXS(float energy) const;
    float getFissionXS(int energy_index) const;
    float getTotalXS(float energy) const;
    float getTotalXS(int energy_index) const;
    float getTransportXS(int energy_index) const;
    float getTransportXS(float energy) const;
    bool usesThermalScattering();
    bool isRescaled() const;
    int getEnergyGridIndex(float energy) const;

    /* IMPORTANT: The following eight class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be transformed
     * into numpy arrays */
    void retrieveXSEnergies(float* energies, int num_xs, 
                                            char* xs_type) const;
    void retrieveXS(float* xs, int num_xs, char* xs_type) const;


	void setElasticXS(double* energies, int num_energies,
                       double* elastic_xs, int num_xs);
	void setCaptureXS(double* energies, int num_energies,
                       double* capture_xs, int num_xs);
	void setFissionXS(double* energies, int num_energies,
                       double* fission_xs, int num_xs);

	void setMultigroupElasticXS(double* energies, int num_energies, 
                                double* elastic_xs, int num_xs);
	void setMultigroupCaptureXS(double* energies, int num_energies, 
                                double* capture_xs, int num_xs);
	void setMultigroupFissionXS(double* energies, int num_energies, 
                                double* fission_xs, int num_xs);


    void setA(int A);
    void setAO(float AO);
    void setN(float N);
    void setTemperature(float T);
	void neglectThermalScattering();
	void useThermalScattering();

    Isotope* clone();

    void sampleCollisionType(neutron* neutron);
    void collideNeutron(neutron* neutron);
    float getDistanceTraveled(neutron *neutron);
    float getThermalScatteringEnergy(float energy);

    int getNumThermalCDFs();
    int getNumThermalCDFBins();
    void retrieveThermalCDFs(float* cdfs, int num_values);
    void retrieveThermalDistributions(float* cdfs, int num_values);
    void retrieveEtokT(float* E_to_kT, int num_cdfs);
    void retrieveEprimeToE(float* Eprime_to_E, int num_bins);
};


/**
 * This method returns the index for a certain energy (eV) into
 * the uniform lethargy grid if this Isotope's
 * cross-sections have been rescaled
 * @param energy the energy (eV) of interest
 * @return the index into the uniform lethargy grid
 */
inline int Isotope::getEnergyGridIndex(float energy) const {

	int index;

	energy = log10(energy);

	if (energy > _end_energy)
		index = _num_energies - 1;
	else if (energy < _start_energy)
		index = 0;
	else
		index = int(floor((energy - _start_energy) / _delta_energy));

	return index;
}

#endif

#endif /* ISOTOPE_H_ */
