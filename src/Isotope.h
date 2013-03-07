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

#include <vector>
#include <map>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include "interpolate.h"
#include "integrate.h"
#include "arraycreator.h"
#include "xsreader.h"
#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include "Tally.h"

/* Types of collisions */
typedef enum collisionTypes{
	ELASTIC,
	ABSORPTION,
	CAPTURE,
	FISSION,
	TRANSPORT,
	DIFFUSION,
	LEAKAGE,
	TOTAL
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

    void setIsotopeType(char* isotope);
    void setA(int A);
    void setAO(float AO);
    void setN(float N);

    void setTemperature(float T);

	void loadXS(char* filename, collisionType type);
	void setElasticXS(float* elastic_xs, float* elastic_xs_energies,
			  int num_elastic_xs, scatterAngleType type);
	void setElasticAngleType(scatterAngleType type);
	void setAbsorptionXS(float* absorb_xs, float* absorb_xs_energies,
			     int num_absorb_xs);
	void setCaptureXS(float* capture_xs, float* capture_xs_energies,
			     int num_capture_xs);
	void setFissionXS(float* fission_xs, float* fission_xs_energies,
			  int num_fission_xs);
	void generateCaptureXS();
	void rescaleXS(float* new_energies, int num_energies);
	Isotope* clone();

	collisionType getCollisionType(float energy);
	collisionType collideNeutron(float energy);

	float getThermalScatteringEnergy(float energy);
	void initializeThermalScattering(float start_energy, float end_energy,
					 int num_bins, int num_distributions);
	float thermalScatteringProb(float E_prime_to_E, int dist_index);

	void addTally(Tally *tally);
	void clearTallies();

};

#endif /* ISOTOPE_H_ */
