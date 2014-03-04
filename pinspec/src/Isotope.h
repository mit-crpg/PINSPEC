/**
 * @file Isotope.h
 * @brief The Isotope class.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 13, 2012
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
#include "vector.h"
#include "Neutron.h"
#endif

/**
 * @class Isotope Isotope.h "pinspec/src/Isotope.h"
 * @brief The Isotope represents a nuclide at some temperature.
 * @details The Isotope class represents a nuclide and all of its properties
 *          which are relevant to reactor physics calculations.
 */
class Isotope {	
private:
    /** The name of the isotope-periodic table name followed by atomic number */
    char* _isotope_name;
    /** A static class variable to generate a UID for each new isotope */
    static int _n;
    /** The isotope's unique identifier */
    int _uid;
    /** Atomic number */
    int _A;
    /** Atomic number squared (an optimization for speedup) */
    int _A_squared;
    /** Atomic number plus one squared (an optimization for speeup):
     *  \f$ (A+1)^2 \f$ */
    int _A_plus_one_squared;
    /** \f$ \alpha = \left(\frac{A-1}{A+1}\right)^2\f$ */
    float _alpha;
    /** \f$ \eta = \left(\frac{A+1}{2\sqrt{A}}\right)^2 \f$ */
    float _eta;
    /** \f$ \rho =  \left(\frac{A+1}{2\sqrt{A}}\right) \f$ */
    float _rho;
    /** Atomic number ratio */
    float _AO;
    /** Temperature of the isotope in degrees Kelvin */
    float _T;
    /** The average cosine of the scattering angle: 
     * \f$ \left<\mu\right> = \frac{2}{3A} \f$ */
    float _mu_avg;
    /** Whether isotope is fissionable or not */
    bool _fissionable;
    /** Whether cross-sections are rescaled on uniform lethargy grid */
    bool _rescaled;

    /** The number of elastic scattering cross-section data points */
    int _num_elastic_xs;
    /** Array of microscopic elastic scattering cross-section values */
    float* _elastic_xs;
    /** Array of elastic scattering cross-section energies (eV) */
    float* _elastic_xs_energies;
    /** Whether or not the elastic scattering cross-section is rescaled onto
     * a uniform lethargy grid */
    bool _elastic_rescaled;

    /** The number of capture cross-section data points */
    int _num_capture_xs;
    /** Array of microscopic capture cross-section values */
    float* _capture_xs;
    /** Array of capture cross-section energies (eV) */
    float* _capture_xs_energies;
    /** Whether or not the capture cross-section is rescaled onto
     * a uniform lethargy grid */
    bool _capture_rescaled;

    /** The number of fission cross-section data points */
    int _num_fission_xs;
    /** Array of microscopic fission cross-section values */
    float* _fission_xs;
    /** Array of fission cross-section energies (eV) */
    float* _fission_xs_energies;
    /** Whether or not the fission cross-section is rescaled onto
     * a uniform lethargy grid */
    bool _fission_rescaled;

    /** The number of absorption cross-section data points */
    int _num_absorb_xs;
    /** Array of microscopic absorption cross-section values */
    float* _absorb_xs;
    /** Array of absorption cross-section energies (eV) */
    float* _absorb_xs_energies;

    /** The number of total cross-section data points */
    int _num_total_xs;
    /** Array of microscopic total cross-section values */
    float* _total_xs;
    /** Whether or not the total cross-section is rescaled onto
     * a uniform lethargy grid */
    float* _total_xs_energies;

    /** Number of rescaled cross-section values on uniform lethargy grid */
    int _num_energies;
    /** Starting lethargy for uniform lethargy grid */
    float _start_lethargy;
    /** Final lethargy for uniform lethargy grid */
    float _end_lethargy;
    /** Space between lethargies in uniform grid */
    float _delta_lethargy;

    /** Whether or not to use thermal scattering */
    bool _use_thermal_scattering;
    /** The high energy cutoff for the thermal scattering treatment (eV) */
    float _thermal_cutoff;
    /** Boltzmann's constant */
    float _kB;
    /** The number of thermal scattering CDFs */
    int _num_thermal_cdfs;
    /** The number of bins per thermal scattering CDFs */
    int _num_thermal_cdf_bins;
    /** The number of thermal scattering PDFs */
    float* _thermal_dist;
    /** 2D array of thermal scattering CDFs */
    float** _thermal_cdfs;
    /** Array of the \f$ \frac{E}{kT} \f$ values for each PDF/CDF */
    float* _E_to_kT;
    /** Array of \f$ \frac{E}{E'} \f$ for each PDF/CDF */
    float* _Eprime_to_E;

    void loadXS();
    void setElasticXS(float* elastic_xs, float* elastic_xs_energies,								   int num_elastic_xs);
    void setCaptureXS(float* capture_xs, float* capture_xs_energies,
						   int num_capture_xs);
    void setFissionXS(float* fission_xs, float* fission_xs_energies,
			                          int num_fission_xs);
    void rescaleXS(float start_energy, float end_energy, int num_energies);
    void generateAbsorptionXS(float start_energy, float end_energy, 
			      int num_energies);
    void generateTotalXS(float start_energy, float end_energy, 
			 int num_energies);

    void initializeThermalScattering(float start_energy, float end_energy,
					 int num_bins, int num_distributions);
    float thermalScatteringProb(float E_prime_to_E, int dist_index);

public:
    Isotope(char *_isotope_name);
    virtual ~Isotope();

    void parseName(char* isotope_name);
    void makeFissionable();

    char* getIsotopeName() const;
    int getUid() const;
    int getA() const;
    float getAlpha() const;
    float getTemperature() const;
    float getMuAverage() const;
    bool isFissionable() const;
    float getThermalScatteringCutoff();

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

    /* IMPORTANT: The following eight class method prototypes must
     *  not be changed without changing Geometry.i to allow for the 
     * data arrays to be transformed into numpy arrays */
    void retrieveXSEnergies(float* energies, int num_xs, char* xs_type) const;
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


    void loadXS(char* xs_type);
    void setA(int A);
    void setTemperature(float T);
    void neglectThermalScattering();
    void setThermalScatteringCutoff(float cutoff_energy);
    void useThermalScattering();

    Isotope* clone();

    void sampleCollisionType(neutron* neutron);
    float getDistanceTraveled(neutron* neutron);
    void collideNeutron(neutron* neutron);
    float getThermalScatteringEnergy(float energy);

    int getNumThermalCDFs();
    int getNumThermalCDFBins();
    void retrieveThermalCDFs(float* cdfs, int num_values);
    void retrieveThermalPDFs(float* pdfs, int num_values);
    void retrieveEtokT(float* E_to_kT, int num_cdfs);
    void retrieveEprimeToE(float* Eprime_to_E, int num_bins);
};


/**
 * @brief This method returns the index for a certain energy (eV) into
 *        the Isotope's uniform lethargy grid.
 * @details The index computed is that of nearest energy less than
 *          or equal to the input energy.
 * @param energy the energy (eV) of interest
 * @return the index into the uniform lethargy grid
 */
inline int Isotope::getEnergyGridIndex(float energy) const {

    int index;

    /* Compute the lethargy */
    float lethargy = log10(energy);

    /* If the energy is outside of the grid, pin the energy to the max/min */
    if (lethargy > _end_lethargy)
        index = _num_energies - 1;
    else if (lethargy < _start_lethargy)
        index = 0;
    /* If the lethargy is within the grid, comput the index */
    else
        index = int(floor((lethargy - _start_lethargy) / _delta_lethargy));

    return index;
}


#endif /* ISOTOPE_H_ */
