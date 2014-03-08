/**
 * @file Fissioner.h
 * @brief The Fissioner class for fission neutron emission.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 15, 2012
 */

#ifndef FISSIONER_H_
#define FISSIONER_H_

#ifdef __cplusplus
#include "log.h"
#include "integrate.h"
#include "interpolate.h"
#include "arraycreator.h"
#include <random>
using namespace std;
#endif


/**
 *
 * @class Fissioner Fissioner.h "pinspec/src/Fissioner.h" 
 * @brief The Fissioner represents the physics of fission neutron emission.
 * @details The Fissioner class contains a cumulative distribution 
 *          function (CDF) for the chi spectrum of fission neutron energies. The
 *          Fissioner can stochastically sample from the CDF to generate 
 *          the neutron fission emisson spectrum. The Watt spectrum given in 
 *          the "Fundamentals of Nuclear Reactor Physics", E. E. Lewis is
 *          used to generate the CDF:
 *          
 *          \f$ \chi = 0.453 * exp(-1.036E) * sinh(\sqrt{2.29E}) \f$
 */
class Fissioner {

private:
    /** The number of Watt spectrum CDF bins */
    int _num_bins;
    /** The array of CDF values */
    float* _cdf;
    /** The array of CDF energies */
    float* _cdf_energies;
    /** The maximum fission emission energy in MeV for the CDF */
    float _E_max;
    /** The random number seed */
    unsigned int _seed;

public:
    Fissioner();
    virtual ~Fissioner();
    int getNumBins();
    void setNumBins(int num_bins);
    void setEMax(float E_max);
    void setRandomNumberSeed(unsigned int seed);
    void initializeRandomNumberGenerator();

    void buildCDF();
    float wattSpectrum(float energy);
    float emitNeutronMeV();
    float emitNeutroneV();
    
    void retrieveCDF(float* cdf, int num_bins);
    void retrieveCDFEnergies(float* cdf_energies, int num_bins);
};


#endif /* FISSIONER_H_ */
