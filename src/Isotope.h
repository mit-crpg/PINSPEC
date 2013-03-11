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
#endif
#include "interpolate.h"
#include "integrate.h"
#include "arraycreator.h"
#include "xsreader.h"
#include "log.h"
#include "Tally.h"
#include "Neutron.h"


/* Types of collisions */
typedef enum collisionTypes{
	ELASTIC,
	CAPTURE,
	FISSION,
	LEAKAGE,
} collisionType;

/* Types of angular scattering distributions */
typedef enum scatterAngleTypes {
	ISOTROPIC_CM,
	ISOTROPIC_LAB
} scatterAngleType;


/**
 * The isotope class represents an isotope and all of its properties
 * which are relevant to neutronics
 */
#ifdef __cplusplus
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

	int _num_elastic_xs;
	float* _elastic_xs;
	float* _elastic_xs_energies;
	scatterAngleType _elastic_angle;
	int _num_absorb_xs;
	float* _absorb_xs;
	float* _absorb_xs_energies;
	int _num_capture_xs;
	float* _capture_xs;
	float* _capture_xs_energies;
	int _num_fission_xs;
	float* _fission_xs;
	float* _fission_xs_energies;
	int _num_total_xs;
	float* _total_xs;
	float* _total_xs_energies;

	/* Values related to rescaled cross-sections on a uniform energy grid */
	bool _rescaled;
	binSpacingTypes _scale_type;
	int _num_energies;
	float _start_energy;
	float _end_energy;
	float _delta_energy;

	/* Map of keys (xs types) with values (getXS functions for xs types) */
	std::map<collisionType, float(Isotope::*)(float) const> _xs_handles;
	std::vector<Tally*> _tallies;

	int _num_thermal_cdfs;
	int _num_thermal_cdf_bins;
	float* _thermal_dist;
	float** _thermal_cdfs;
	float* _E_to_kT;
	float* _Eprime_to_E;
	float _kB;

	void loadXS();
	void setElasticXS(float* elastic_xs, float* elastic_xs_energies,
			  int num_elastic_xs, scatterAngleType type);
	void setElasticAngleType(scatterAngleType type);
	void setCaptureXS(float* capture_xs, float* capture_xs_energies,
			     int num_capture_xs);
	void setFissionXS(float* fission_xs, float* fission_xs_energies,
			  int num_fission_xs);
	void rescaleCrossSections(float start_energy, float end_energy,
								int num_energies, binSpacingTypes scale_type);	
	void rescaleXS(float* new_energies, int num_energies);


	void initializeThermalScattering(float start_energy, float end_energy,
					 int num_bins, int num_distributions);
	float thermalScatteringProb(float E_prime_to_E, int dist_index);

	void clearTallies();

public:
	Isotope(char *_isotope_name);
    virtual ~Isotope();

	void parseName();
	void makeFissionable();

	char* getIsotopeType() const;
    int getA() const;
    float getAlpha() const;
    float getN() const;
    float getAO() const;
    float getTemperature() const;
    float getMuAverage() const;
	bool isFissionable() const;

    int getNumXSEnergies() const;
    float getElasticXS(float energy) const;
    float getElasticXS(int energy_index) const;
    scatterAngleType getElasticAngleType() const;
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
	bool isRescaled();
	int getEnergyGridIndex(float energy);

    /* IMPORTANT: The following two class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be transformed
     * into numpy arrays */
    void retrieveXSEnergies(float* energies, int num_xs) const;
    void retrieveXS(float* xs, int num_xs, char* xs_type) const;

    void setIsotopeType(char* isotope);
    void setA(int A);
    void setAO(float AO);
    void setN(float N);
    void setTemperature(float T);
    void setNumBatches(int num_batches);

	Isotope* clone();

	collisionType getCollisionType(float energy);
	collisionType collideNeutron(neutron* neutron);
	float getDistanceTraveled(neutron *neutron);
	float getThermalScatteringEnergy(float energy);

    int getNumThermalCDFs();
    int getNumThermalCDFBins();
    void retrieveThermalCDFs(float* cdfs, int num_values);
    void retrieveThermalDistributions(float* cdfs, int num_values);
    void retrieveEtokT(float* E_to_kT, int num_cdfs);
    void retrieveEprimeToE(float* Eprime_to_E, int num_bins);
    
	void addTally(Tally *tally);

    void computeBatchStatistics();
    void computeScaledBatchStatistics(float scale_factor);
    void outputBatchStatistics(char* directory, char* suffix);
};

#endif

#endif /* ISOTOPE_H_ */
