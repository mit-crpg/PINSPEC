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

#include <map>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include "interpolate.h"
#include "integrate.h"
#include "arraycreator.h"
#include "xsreader.h"
#include "log.h"

/* Types of collisions */
typedef enum collisionTypes{
	CAPTURE,
	ELASTIC,
	FISSION,
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
class Isotope
{
private:
	char* _isotope_name;
	int _A;
	float _alpha;
	float _eta;
	float _rho;
	float _N;
	float _T;
	float _mu_avg;
	float* _capture_xs;
	float* _capture_xs_energies;
	int _num_capture_xs;
	float* _scatter_xs;
	float* _scatter_xs_energies;
	int _num_scatter_xs;
	scatterAngleType _scatter_angle;
	float* _elastic_xs;
	float* _elastic_xs_energies;
	int _num_elastic_xs;
	scatterAngleType _elastic_angle;
	float* _fission_xs;
	float* _fission_xs_energies;
	int _num_fission_xs;
	float* _total_xs;
	float* _total_xs_energies;
	int _num_total_xs;

	/* Map of keys (xs types) with values (getXS functions for xs types) */
	std::map<collisionType, float(Isotope::*)(float) const> _xs_handles;

	int _num_thermal_cdfs;
	int _num_thermal_cdf_bins;
	float* _thermal_dist;
	float** _thermal_cdfs;
	float* _E_to_kT;
	float* _Eprime_to_E;
	float _kB;
public:
	Isotope();
    virtual ~Isotope();

    char* getIsotopeType() const;
    int getA() const;
    float getAlpha() const;
    float getN() const;
    float getTemperature() const;
    float getMuAverage() const;
    float getCaptureXS(float energy) const;
    float getCaptureXS(int energy_index) const;
    float getElasticXS(float energy) const;
    float getElasticXS(int energy_index) const;
    scatterAngleType getElasticAngleType() const;
    float getFissionXS(float energy) const;
    float getFissionXS(int energy_index) const;
    float getAbsorbXS(float energy) const;
    float getAbsorbXS(int energy_index) const;
    float getTotalXS(float energy) const;
    float getTotalXS(int energy_index) const;
    float getTransportXS(int energy_index) const;
    float getTransportXS(float energy) const;
    bool usesThermalScattering();

    void setIsotopeType(char* isotope);
    void setA(int A);
    void setN(float N);
    void setTemperature(float T);
    void loadXS(char* filename, collisionType type, char* delimiter);
    void setCaptureXS(float* capture_xs, float* capture_xs_energies,
    											int num_capture_xs);
    void setElasticXS(float* elastic_xs, float* elastic_xs_energies,
							int num_elastic_xs, scatterAngleType type);
    void setElasticAngleType(scatterAngleType type);
    void setFissionXS(float* fission_xs, float* fission_xs_energies,
												int num_fission_xs);

    collisionType getCollisionType(float energy);
    float getThermalScatteringEnergy(float energy);
    Isotope* clone();
    void initializeThermalScattering(float start_energy, float end_energy,
    								int num_bins, int num_distributions);
    float thermalScatteringProb(float E_prime_to_E, int dist_index);
    void rescaleXS(float* new_energies, int num_energies);
};

#endif /* ISOTOPE_H_ */
