#include "Fissioner.h"


/**
 * @brief Default Fissioner constructor.
 * @details Initializes the Fissioner with 1E5 bins and a maximum emission 
 *          energy of 20 MeV.
 */
Fissioner::Fissioner() {
    _num_bins = 100000;
    _E_max = 20;
    buildCDF();
}


/**
 * @brief Fissioner destructor deallocates memory for CDF arrays.
 */
Fissioner::~Fissioner() {
    if (_num_bins != 0) {
         delete [] _cdf;
         delete [] _cdf_energies;
    }
}


/**
 * @brief Returns the number of bins that the Fissioner uses for the CDF.
 * @return the number of CDF bins
 */
int Fissioner::getNumBins() {
    return _num_bins;
}


/**
 * @brief Sets the number of bins that we wish to use for the CDF.
 * @param num_bins the number of CDF bins
 */
void Fissioner::setNumBins(int num_bins) {
    _num_bins = num_bins;
}


/**
 * @brief Sets the maximum energy value for the CDF.
 * @param E_max the maximum CDF energy value in MeV
 */
void Fissioner::setEMax(float E_max) {
    _E_max = E_max;
}


/**
 * @brief Sets the random number seed for collision probability sampling.
 * @details This method is called by the Geometry::setRandomNumberSeed(...)
 *          method and should not be called directly by the user.
 * @param seed the random number seed
 */
void Fissioner::setRandomNumberSeed(unsigned int seed) {
    _seed = seed;
}


/**
 * @brief Initializes the random number seed for random number sampling.
 * @details This method is called by the 
 *          Geometry::initializeRandomNumberGenerator()
 *          method and should not be called directly by the user. This
 *          method likewise calls the 
 *          Material::initializeRandomNumberGenerator()
 *          method for its material.
 */
void Fissioner::initializeRandomNumberGenerator() {
    srand(_seed);

    log_printf(NORMAL, "Initializing fissioner's random number seed to %d", _seed);

    log_printf(NORMAL, "First random #: %f\n", float(rand()) / RAND_MAX);


    default_random_engine generator(_seed);
    uniform_real_distribution<double> distribution(0.0,100., generator());

    log_printf(NORMAL, "my first random number %f\n", distribution(generator));
}


/**
 * @brief Builds the CDF of the Watt Spectrum.
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
 * @brief Returns the chi value for a given energy from the Watt spectrum.
 * @param energy a fission energy value in MeV
 * @return the value of chi \f$\chi\f$ at that energy
 */
float Fissioner::wattSpectrum(float energy) {
    return (0.453 * exp(-1.036 * energy) * sinh(sqrt(2.29 * energy)));
}


/**
 * @brief Sample the Fissioner's CDF and return a neutron energy from fission
 * @return a neutron energy in MeV
 * @see Fissioner::emitNeutroneV()
 */
float Fissioner::emitNeutronMeV() {

    /* Check that the CDF has been built */
    if (_num_bins == 0)
        log_printf(ERROR, "Unable to sample Fissioner CDF since it "
	 			"has not yet been created");

    /* Return an interpolate value from the fission energy CDF */
    return linearInterp<float, float, float>(_cdf, _cdf_energies, _num_bins,
					     float(rand()) / RAND_MAX);
}


/**
 * @brief Sample the Fissioner's CDF and return a neutron energy from fission
 * @return a neutron energy in eV
 * @see Fissioner::emitNeutronMeV()
 */
float Fissioner::emitNeutroneV() {
    return emitNeutronMeV() * 1E6;
}


/**
 * @brief Fills an input array with the Fissioner's CDF values
 * @details This class method is a helper function to allow the user 
 *          to access the CDF array through in Python via SWIG by inputting 
 *          a numpy array to this function. For example, a user may use this 
 *          function in Python as follows:
 * 
 * @code
 *     num_cdf_bins = fissioner.getNumBins()
 *     numpy_cdf = numpy.zeros(num_cdf_bins)
 *     fissioner.retrieveCDF(numpy_cdf)
 * @endcode
 *     
 * @param cdf an array of the same length as the Fissioner's CDF
 * @param num_bins the length of the input array for the CDF
 */
void Fissioner::retrieveCDF(float* cdf, int num_bins) {

    /* Check that the input array is the correct size */
    if (num_bins != _num_bins)
        log_printf(ERROR, "Unable to retrieve the Fissioner's CDF into an"
                        " array of length %d since the CDF has %d bins",
                                                    num_bins, _num_bins);    

    /* Fill the input array with the CDF */
    for (int i=0; i < _num_bins; i++)
        cdf[i] = _cdf[i];        
}



/**
 * @brief Fills an input array with the Fissioner's CDF energies
 * @details This class method is a helper function to allow the user to access 
 *          the CDF array through in Python via SWIG by inputting a numpy array 
 *          to this function. For example, a user may use this function in 
 *          Python as follows:
 * 
 * @code
 *     num_cdf_energies = fissioner.getNumBins()
 *     numpy_cdf_energies = numpy.zeros(num_cdf_energies)
 *     fissioner.retrieveCDF(numpy_cdf_energies)
 * @endcode
 *     
 * @param cdf_energies an array of the same length as the Fissioner's CDF
 * @param num_bins the length of the input array for the CDF
 */
void Fissioner::retrieveCDFEnergies(float* cdf_energies, int num_bins) {

    /* Check that the input array is the correct size */
    if (num_bins != _num_bins)
        log_printf(ERROR, "Unable to retrieve the Fissioner's CDF energies "
                        "into an array of length %d since the CDF has %d bins",
                                                        num_bins, _num_bins);    
    /* Fill the input array with the CDF energies */    
    for (int i=0; i < _num_bins; i++)
       cdf_energies[i] = _cdf_energies[i]; 
}
