#include "Isotope.h"


int Isotope::_n = 1;


/**
 * @brief Isotope constructor.
 * @details Searches the cross-section library for appropriately named
 *          files using the isotope name and loads the capture, scatter, and 
 *          fission (if file is found) cross-sectoins. By default, the
 *          constructor rescales the cross-sections onto a uniform lethargy
 *          grid of 100,000 values between 1E-5 eV and 20 MeV. In addition,
 *          the constructor creates thermal scattering CDFs for the isotope
 *          at 300K by default.
 */
Isotope::Isotope(char* isotope_name){

    /* Copies the isotope's name and uses it to find the cross-section files */
    parseName(isotope_name);

    _uid = _n;
    _n++;

    _A_squared = _A * _A;
    _A_plus_one_squared = (_A + 1) * (_A + 1);

    _T = 300;
    _kB = 8.617332E-5;             /* boltzmann's constant (ev / K) */
    _fissionable = false;

    /** Set the default random number seed */
    setRandomNumberSeed(SEED);

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
    _start_lethargy = 1E-5;
    _end_lethargy = 2E7;
    _num_energies = 100000;
    rescaleXS(_start_lethargy, _end_lethargy, _num_energies);
 
    /* By default the thermal scattering cdfs have not been initialized */
    _thermal_cutoff = 4.0;
    _use_thermal_scattering = true;
    _num_thermal_cdfs = 0;
    _num_thermal_cdf_bins = 0;

    /* FIXME: may need to update these defaults since they may need to
     * depend on the isotope's atomic number */
    initializeThermalScattering(1E-6, 15, 1000, 15);
}



/**
 * @brief Isotope destructor deletes arrays of cross-section values that
 *        have been assigned to this isotope.
 */
Isotope::~Isotope() {

    delete [] _isotope_name;

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
 * @brief Parse input name and set atomic number for isotope.
 */
void Isotope::parseName(char* isotope_name){

    int A = 0;
    int i = 0;
    int length = strlen(isotope_name);

    /* Determine the length of the Isotope's name */
    while (isotope_name[i] != '\0'){

        /* if - encountered, get the A value */
        if (isotope_name[i] == '-') {
            i++;
            A = atoi(&isotope_name[i]);
	    break;
        }

        i++;
    }

    /* Initialize a character array for the Isotope's name */
    _isotope_name = new char[length+1];

    /* Copy the input character array isotope name to the class attribute name */
    for (int k=0; k <= length; k++)
        _isotope_name[k] = isotope_name[k];

    if (A < 1 || A > 300)
	log_printf(ERROR, "Isotope identifier %s is not formatted correctly", 
		   isotope_name);

    /* Set the atomic number of isotope */
    setA(A);

    log_printf(DEBUG, "Isotope %s has atomic number %i", _isotope_name, _A);
}



/**
 * @brief Returns the name of the of isotope.
 * @return character array with name of isotope
 */
char* Isotope::getIsotopeName() const {
    return _isotope_name;
}


/**
 * @brief Returns the unique ID auto-generated for the isotope.
 * @return a unique ID for the isotope
 */
int Isotope::getUid() const {
    return _uid;
}


/**
 * @brief Returns the atomic number of this isotope.
 * @return the atomic number
 */
int Isotope::getA() const {
    return _A;
}


/**
 * @brief Returns the alpha \f$ \alpha = ((A-1)/(A+1))^2 \f$ values for 
 *        this isotope.
 * @return alpha
 */
float Isotope::getAlpha() const {
    return _alpha;
}


/**
 * @brief Return the temperature in degrees Kelvin for this isotope.
 * @return the temperature of this isotope
 */
float Isotope::getTemperature() const {
    return _T;
}


/**
 * @brief Return the average value of the cosine of the scattering angle.
 * @details Returns the average value of cosine of theta for this isotope
 *        in a scattering collision: \f$ \left<\mu\right> = \frac{2}{3A} \f$
 * @return the average for mu
 */
float Isotope::getMuAverage() const {
    return _mu_avg;
}


/**
 * @brief Return whether this isotope is fissionable (true) or not (false)
 * @return whether this isotope is fissionable
 */
bool Isotope::isFissionable() const {
    return _fissionable;
}


/**
 * @brief Return the thermal scattering high energy cutoff.
 * @return the thermal scattering high energy cutoff (eV)
 *
 */
float Isotope::getThermalScatteringCutoff() {
    return _thermal_cutoff;
}


/**
 * @brief Returns the number of energies for a particular cross-section type
 * @details Returns the number of energies for 'capture', 'elastic', 'fission',
 *          'fission' and 'absorption cross-section types.
 * @param xs_type a character array for the xs_type
 * @return 
 */
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


/**
 * @brief Fills an array with cross-section energy values.
 * @details This method is a helper function to allow PINSPEC users to
 *          get access to the isotope's nuclear data in Python. A user
 *          must initialize a numpy array of the correct size (ie, 
 *          a float64 array the length of the number of cross-section
 *          values) as input to this function. This function then fills
 *          the numpy array with the energy values for the isotope's
 *          cross-section data. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          num_energies = isotope.getNumXSEnergies()
 *          energies = numpy.zeros(num_energies)          
 *          isotope.retrieveXSEnergies(energies, num_energies, 'capture')
 * @endcode
 * 
 * @param energies an array to fill with the cross-section energies
 * @param num_xs the number of cross-section values
 * @param xs_type the type of cross-section
 */
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


/**
 * @brief Fills an array with microscopic cross-section values.
 * @details This method is a helper function to allow PINSPEC users to
 *          get access to the isotope's nuclear data in Python. A user
 *          must initialize a numpy array of the correct size (ie, 
 *          a float64 array the length of the number of cross-section
 *          values) as input to this function. This function then fills
 *          the numpy array with the data values for one of the isotope's
 *          cross-sections. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          num_xs = isotope.getNumXSEnergies()
 *          xs = numpy.zeros(num_xs)          
 *          isotope.retrieveXS(xs, num_xs, 'capture')
 * @endcode
 * 
 * @param xs an array to fill with the microscopic cross-section data
 * @param num_xs the number of cross-section values
 * @param xs_type the type of cross-section
 */
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


/**
 * @brief Sets the elastic cross-section data for this isotope.
 * @details This is a helper function for users to assign the cross-section
 *          data for an isotope from a numpy array in Python. Although the 
 *          prototype for this function seems to require four arguments - 
 *          two with arrays of data for energies and cross-sections and two
 *          for the length of each array - in Python one must only give the
 *          method a handle to each of two arrays. A user may call this 
 *          method from within Python as follows:
 *
 * @code
 *          energies = numpy.array([1E-3, 0.1, 1., 10., 1000.])
 *          xs = numpy.array([1000., 1000., 10., 1., 0.1])
 *          isotope.setElasticXS(energies, xs)
 * @endcode
 * 
 * @param energies an array of energy values
 * @param num_energies the number of data points
 * @param elastic_xs the microscopic elastic cross-section values
 * @param num_xs the number of cross-section values
 */
void Isotope::setElasticXS(double* energies, int num_energies,
                               double* elastic_xs, int num_xs) {

    /* Check that the # of xs values is 1 less than the # energy bounds */
    if (num_xs != num_energies)
		log_printf(ERROR, "Unable to set elastic xs for isotope %s "
			   "since the number of xs values is %d while the "
			   "number of energies is %d", _isotope_name, num_xs, 
			   num_energies);

    /* Check that all energies are >0 and monotonically increasing */
    double prev_energy = energies[0];

    for (int i=0 ; i < num_xs; i++) {

	/* Check for monotonically increasing energy */
	if (energies[i] < prev_energy)
		log_printf(ERROR, "Unable to set elastic xs for isotope"
			   " %s since all xs energies must be monotonically "
					"increasing", _isotope_name);
		/* Check that energy is non-negative */		
	if (energies[i] < 0.0)
		log_printf(ERROR, "Unable to set elastic xs for isotope" 
			"%s since all xs energies must be "
			"non-negative", _isotope_name);
    }

	/* Check that all xs are non-negative */
    for (int i=0; i < num_xs; i++) {
        if (elastic_xs[i] < 0.0)
	    log_printf(ERROR, "Unable to set elastic xs for isotope" 
		       "%s since all xs values must be "
			"non-negative", _isotope_name);
    }

    log_printf(INFO, "Setting elastic xs for isotope %s", _isotope_name);
	
    /* Free memory for old elastic xs from ENDF */
    delete [] _elastic_xs;
    delete [] _elastic_xs_energies;

    /* Allocate memory for the new elastic xs and energies */
    _num_elastic_xs = num_xs;
    _elastic_xs = new float[_num_elastic_xs];
    _elastic_xs_energies = new float[_num_elastic_xs];

    for (int i=0; i < _num_elastic_xs; i++) {
        _elastic_xs[i] = elastic_xs[i];
        _elastic_xs_energies[i] = energies[i];
    }

    _elastic_rescaled = false;

    rescaleXS(pow(10., _start_lethargy), pow(10., _end_lethargy), _num_energies);

    return;
}



/**
 * @brief Sets the capture cross-section data for this isotope.
 * @details This is a helper function for users to assign the cross-section
 *          data for an isotope from a numpy array in Python. Although the 
 *          prototype for this function seems to require four arguments - 
 *          two with arrays of data for energies and cross-sections and two
 *          for the length of each array - in Python one must only give the
 *          method a handle to each of two arrays. A user may call this 
 *          method from within Python as follows:
 *
 * @code
 *          energies = numpy.array([1E-3, 0.1, 1., 10., 1000.])
 *          xs = numpy.array([1000., 1000., 10., 1., 0.1])
 *          isotope.setCaptureXS(energies, xs)
 * @endcode 
 *
 * @param energies an array of energy values
 * @param num_energies the number of data points
 * @param capture_xs the microscopic elastic cross-section values
 * @param num_xs the number of cross-section values
 */
void Isotope::setCaptureXS(double* energies, int num_energies,
                               double* capture_xs, int num_xs) {

    /* Check that the # of xs values is 1 less than the # energy bounds */
    if (num_xs != num_energies)
         log_printf(ERROR, "Unable to set capture xs for isotope %s since"
		    " the number of xs values is %d while the number of "
		    "energies is %d", _isotope_name, num_xs, num_energies);

    /* Check that all energies are >0 and monotonically increasing */
    double prev_energy = energies[0];

    for (int i=0 ; i < num_xs; i++) {

        /* Check for monotonically increasing energy */
	if (energies[i] < prev_energy)
		log_printf(ERROR, "Unable to set capture xs for isotope %s"
                        " since all xs energies must be monotonically "
					"increasing", _isotope_name);
	/* Check that energy is non-negative */		
	if (energies[i] < 0.0)
		log_printf(ERROR, "Unable to set capture xs for isotope %s"
			   " since all xs energies must be "
			   "non-negative", _isotope_name);
	}

    /* Check that all xs are non-negative */
    for (int i=0; i < num_xs; i++) {

        if (capture_xs[i] < 0.0)
	    log_printf(ERROR, "Unable to set capture xs for isotope %s"
		       " since all xs values must be "
		       "non-negative", _isotope_name);
	}

    log_printf(INFO, "Setting capture xs for isotope %s", _isotope_name);
	
    /* Free memory for old capture xs from ENDF */
    delete [] _capture_xs;
    delete [] _capture_xs_energies;
  
    /* Allocate memory for the new capture xs and energies */
    _num_capture_xs = num_xs;
    _capture_xs = new float[_num_capture_xs];
    _capture_xs_energies = new float[_num_capture_xs];

    for (int i=0; i < _num_capture_xs; i++) {
        _capture_xs[i] = capture_xs[i];
        _capture_xs_energies[i] = energies[i];
    }

    _capture_rescaled = false;

    rescaleXS(pow(10., _start_lethargy), pow(10., _end_lethargy), 
	      _num_energies);

    return;
}



/**
 * @brief Sets the fission cross-section data for this isotope.
 * @details This is a helper function for users to assign the cross-section
 *          data for an isotope from a numpy array in Python. Although the 
 *          prototype for this function seems to require four arguments - 
 *          two with arrays of data for energies and cross-sections and two
 *          for the length of each array - in Python one must only give the
 *          method a handle to each of two arrays. A user may call this 
 *          method from within Python as follows:
 *
 * @code
 *          energies = numpy.array([1E-3, 0.1, 1., 10., 1000.])
 *          xs = numpy.array([1000., 1000., 10., 1., 0.1])
 *          isotope.setFissionXS(energies, xs)
 * @endcode 
 *
 * @param energies an array of energy values
 * @param num_energies the number of data points
 * @param fission_xs the microscopic fission cross-section values
 * @param num_xs the number of cross-section values
 */
void Isotope::setFissionXS(double* energies, int num_energies,
                               double* fission_xs, int num_xs) {

    /* Check that the # of xs values is 1 less than the # energy bounds */
    if (num_xs != num_energies)
	log_printf(ERROR, "Unable to set fission xs for isotope %s since"
		   " the number of xs values is %d while the number of "
		   "energies is %d", _isotope_name, num_xs, num_energies);

    /* Check that all energies are >0 and monotonically increasing */
    double prev_energy = energies[0];

    for (int i=0 ; i < num_xs; i++) {

        /* Check for monotonically increasing energy */
	if (energies[i] < prev_energy)
		log_printf(ERROR, "Unable to set fission xs for isotope %s"
			   " since all xs energies must be monotonically "
			   "increasing", _isotope_name);
	/* Check that energy is non-negative */		
	if (energies[i] < 0.0)
		log_printf(ERROR, "Unable to set fission xs for isotope %s"
			   " since all xs energies must be "
			   "non-negative", _isotope_name);
	}

    /* Check that all xs are non-negative */
    for (int i=0; i < num_xs; i++) {
	if (fission_xs[i] < 0.0)
		log_printf(ERROR, "Unable to set fission xs for isotope %s"
			   " since all xs values must be "
			   "non-negative", _isotope_name);
    }

    log_printf(INFO, "Setting fission xs for isotope %s", _isotope_name);
	
    /* Free memory for old elastic xs from ENDF */
    delete [] _fission_xs;
    delete [] _fission_xs_energies;

    /* Allocate memory for the new fission xs and energies */
    _num_fission_xs = num_xs;
    _fission_xs = new float[_num_fission_xs];
    _fission_xs_energies = new float[_num_fission_xs];

    for (int i=0; i < _num_fission_xs; i++) {
        _fission_xs[i] = fission_xs[i];
        _fission_xs_energies[i] = energies[i];
    }

    _fission_rescaled = false;

    rescaleXS(pow(10., _start_lethargy), pow(10., _end_lethargy), 
	      _num_energies);

    return;
}


/**
 * @brief Sets the multigroup elastic cross-section data for this isotope.
 * @details This is a helper function for users to assign the cross-section
 *          data for an isotope from a numpy array in Python. Although the 
 *          prototype for this function seems to require four arguments - 
 *          two with arrays of data for energies and cross-sections and two
 *          for the length of each array - in Python one must only give the
 *          method a handle to each of two arrays. A user may call this 
 *          method from within Python as follows:
 *
 * @code
 *          energies = numpy.array([1E-3, 0.1, 1., 10., 1000., 1E7])
 *          xs = numpy.array([1000., 1000., 10., 1., 0.1])
 *          isotope.setMultigroupElasticXS(energies, xs)
 * @endcode 
 *
 * @param energies an array of energy bounds
 * @param num_energies the number of energy bounds
 * @param elastic_xs the microscopic elastic multigroup cross-sections
 * @param num_xs the number of multigroup cross-sections
 */
void Isotope::setMultigroupElasticXS(double* energies, int num_energies, 
                                        double* elastic_xs, int num_xs) {

    /* Check that the # of xs values is 1 less than the # energy bounds */
    if (num_xs != num_energies-1)
        log_printf(ERROR, "Unable to set multigroup elastic xs for "
		   "isotope %s since the number of xs values is %d while "
		   "the number of energies is %d", 
		   _isotope_name, num_xs, num_energies);

        /* Check that all energies are >0 and monotonically increasing */
        double prev_energy = energies[0];

        for (int i=0 ; i < num_xs+1; i++) {

	    /* Check for monotonically increasing energy */
	    if (energies[i] < prev_energy)
	        log_printf(ERROR, "Unable to set multigroup elastic xs for "
			"isotope %s since all xs energies must be "
			"monotonically increasing", _isotope_name);
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

        /* The number of elastic xs values for the multigroup case.
         * We approximate the multigroup xs as a continuous energy xs
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

	/* Rescales the cross-section to a uniform lethargy grid */
        _elastic_rescaled = false;
        rescaleXS(pow(10., _start_lethargy), pow(10., _end_lethargy), 
		  _num_energies);

        return;
}



/**
 * @brief Sets the multigroup capture cross-section data for this isotope.
 * @details This is a helper function for users to assign the cross-section
 *          data for an isotope from a numpy array in Python. Although the 
 *          prototype for this function seems to require four arguments - 
 *          two with arrays of data for energies and cross-sections and two
 *          for the length of each array - in Python one must only give the
 *          method a handle to each of two arrays. A user may call this 
 *          method from within Python as follows:
 *
 * @code
 *          energies = numpy.array([1E-3, 0.1, 1., 10., 1000., 1E7])
 *          xs = numpy.array([1000., 1000., 10., 1., 0.1])
 *          isotope.setMultigroupCaptureXS(energies, xs)
 * @endcode 
 *
 * @param energies an array of energy bounds
 * @param num_energies the number of energy bounds
 * @param capture_xs the microscopic capture multigroup cross-sections
 * @param num_xs the number of multigroup cross-sections
 */
void Isotope::setMultigroupCaptureXS(double* energies, int num_energies,
                                        double* capture_xs, int num_xs) {

    /* Check that the # of xs values is 1 less than the # energy bounds */
    if (num_xs != num_energies-1)
        log_printf(ERROR, "Unable to set multigroup capture xs for "
		   "isotope %s since the number of xs values is %d while "
		   "the number of energies is %d", 
		   _isotope_name, num_xs, num_energies);

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

    /* The number of capture xs values for the multigroup case. 
     * We approximate the multigroup xs as a continuous energy xs
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
        _capture_xs[i] = capture_xs[i/2];
	_capture_xs[i+1] = capture_xs[i/2];
	_capture_xs_energies[i] = energies[i/2] - 1E-3;
	_capture_xs_energies[i+1] = energies[i/2] + 1E-3;
    }

    /* Highest energy xs value */
    _capture_xs_energies[_num_capture_xs-1] = energies[num_xs];
    _capture_xs[_num_capture_xs-1] = capture_xs[num_xs-1];

    /* Rescale the cross-section to a uniform lethargy grid */
    _capture_rescaled = false;
    rescaleXS(pow(10.,_start_lethargy), pow(10., _end_lethargy), _num_energies);
    
    return;
}




/**
 * @brief Sets the multigroup fission cross-section data for this isotope.
 * @details This is a helper function for users to assign the cross-section
 *          data for an isotope from a numpy array in Python. Although the 
 *          prototype for this function seems to require four arguments - 
 *          two with arrays of data for energies and cross-sections and two
 *          for the length of each array - in Python one must only give the
 *          method a handle to each of two arrays. A user may call this 
 *          method from within Python as follows:
 *
 * @code
 *          energies = numpy.array([1E-3, 0.1, 1., 10., 1000., 1E7])
 *          xs = numpy.array([1000., 1000., 10., 1., 0.1])
 *          isotope.setMultigroupFissionXS(energies, xs)
 * @endcode 
 *
 * @param energies an array of energy bounds
 * @param num_energies the number of energy bounds
 * @param fission_xs the microscopic fission multigroup cross-sections
 * @param num_xs the number of multigroup cross-sections
 */
void Isotope::setMultigroupFissionXS(double* energies, int num_energies,
                                     double* fission_xs, int num_xs) {

    /* Check that the # of xs values is 1 less than the # energy bounds */
    if (num_xs != num_energies-1)
        log_printf(ERROR, "Unable to set multigroup fission xs for "
		   "isotope %s since the number of xs values is %d while "
		   "the number of energies is %d", 
		   _isotope_name, num_xs, num_energies);

    /* Check that all energies are >0 and monotonically increasing */
    double prev_energy = energies[0];

    for (int i=0 ; i < num_xs+1; i++) {

        /* Check for monotonically increasing energy */
        if (energies[i] < prev_energy)
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

    /* The number of fission xs values for the multigroup case *
     * We approximate the multigroup xs as a continuous energy xs
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

    /* Rescale the cross-section onto a uniform lethargy grid */
    _fission_rescaled = false;
    rescaleXS(pow(10., _start_lethargy), pow(10., _end_lethargy), 
	      _num_energies);

    return;
}


/**
 * @brief Returns the microscopic elastic scattering cross-section value for 
 *        some energy.
 * @details Uses linear interpolation to compute the cross-section at a
 *          a certain energy (eV).
 * @param energy the energy (eV) of interest
 * @return the microscopic elastic scattering cross-section
 */
float Isotope::getElasticXS(float energy) const {

    if (_num_elastic_xs != 0) {
        /* Use linear interpolation within a uniform lethargy grid */
        if (_elastic_rescaled) {
	    int lower_index = getEnergyGridIndex(energy);
	    double lower_xs = getElasticXS(lower_index);
	    double upper_xs = getElasticXS(lower_index + 1);
	    double delta_xs = upper_xs - lower_xs;
	    double slope = delta_xs / _delta_lethargy;
	    double lower_e = _start_lethargy + _delta_lethargy * lower_index;
	    double xs =  lower_xs + slope * (log10(energy) - lower_e);
	    log_printf(DEBUG, "low = %f, high = %f, xs = %f", lower_xs, 
			   upper_xs, xs);
	    return xs;
        }
	/* Use linear interpolation without a uniform lethargy grid */
	else {
	    return linearInterp<float, float, float>(_elastic_xs_energies,
				   _elastic_xs, _num_elastic_xs, energy);
	}
    }
    /* If this isotope does not have an elastic xs, return 0 */
    else
        return 0.0;
}


/**
 * @brief Returns the microscopic elastic cross-section value for a certain 
 *        energy index in this isotope's rescaled uniform lethargy grid of 
 *        cross-section data.
 * @param energy_index the index into the energy array
 * @return the microscopic elastic cross-section
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
 * @brief Returns the microscopic absorption cross-section value for some energy.
 * @details Uses linear interpolation to compute the cross-section at a
 *          a certain energy (eV).
 * @param energy the energy (eV) of interest
 * @return the microscopic absorption cross-section
 */
float Isotope::getAbsorptionXS(float energy) const {

    if (_num_absorb_xs != 0) {
        /* Use linear interpolation within a uniform lethargy grid */
        if (_rescaled) {
	    int lower_index = getEnergyGridIndex(energy);
	    double lower_xs = getAbsorptionXS(lower_index);
	    double upper_xs = getAbsorptionXS(lower_index + 1);
	    double delta_xs = upper_xs - lower_xs;
	    double slope = delta_xs / _delta_lethargy;
	    double lower_e = _start_lethargy + _delta_lethargy * lower_index;
	    double xs =  lower_xs + slope * (log10(energy) - lower_e);
	    log_printf(DEBUG, "low = %f, high = %f, xs = %f", lower_xs, 
			   upper_xs, xs);
	    return xs;
	}
	/* Use linear interpolation without a uniform lethargy grid */
	else {
	    return linearInterp<float, float, float>(_absorb_xs_energies,
					 _absorb_xs, _num_absorb_xs, energy);
	}
    }
    /* If this isotope does not have an absorption xs, return 0 */
    else
	return 0.0;
}


/**
 * @brief Returns the microscopic absorption cross-section value for a certain 
 *        energy index in this isotope's rescaled uniform lethargy grid of 
 *        cross-section data.
 * @param energy_index the index into the energy array
 * @return the microscopic absorption cross-section
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
 * @brief Returns the microscopic capture cross-section value for some energy.
 * @details Uses linear interpolation to compute the cross-section at a
 *          a certain energy (eV).
 * @param energy the energy (eV) of interest
 * @return the microscopic absorption cross-section
 */
float Isotope::getCaptureXS(float energy) const{

    if (_num_capture_xs != 0) {
        /* Use linear interpolation within a uniform lethargy grid */
        if (_capture_rescaled) {
	    int lower_index = getEnergyGridIndex(energy);
	    double lower_xs = getCaptureXS(lower_index);
	    double upper_xs = getCaptureXS(lower_index + 1);
	    double delta_xs = upper_xs - lower_xs;
	    double slope = delta_xs / _delta_lethargy;
	    double lower_e = _start_lethargy + _delta_lethargy * lower_index;
	    double xs =  lower_xs + slope * (log10(energy) - lower_e);
	    log_printf(DEBUG, "low = %f, high = %f, xs = %f", lower_xs, 
			   upper_xs, xs);
	    return xs;
	}
	/* Use linear interpolation without a uniform lethargy grid */
	else {
	    return linearInterp<float, float, float>(_capture_xs_energies,
			     _capture_xs, _num_capture_xs, energy);
	}
    }
    /* If this isotope does not have a capture xs, return 0 */
    else
        return 0.0;
}


/**
 * @brief Returns the microscopic capture cross-section value for a certain 
 *        energy index in this isotope's rescaled uniform lethargy grid of 
 *        cross-section data.
 * @param energy_index the index into the energy array
 * @return the microscopic capture cross-section
 */
float Isotope::getCaptureXS(int energy_index) const {

    if (energy_index > _num_capture_xs)
	log_printf(ERROR, "Unable to retrieve capture xs for"
		   " isotope %s since the energy index %d is out of"
		   " bounds", _isotope_name, energy_index);

    return _capture_xs[energy_index];
}


/**
 * @brief Returns the microscopic fission cross-section value for some energy.
 * @details Uses linear interpolation to compute the cross-section at a
 *          a certain energy (eV).
 * @param energy the energy (eV) of interest
 * @return the microscopic fission cross-section
 */
float Isotope::getFissionXS(float energy) const{

    if (_num_fission_xs != 0) {
        /* Use linear interpolation within a uniform lethargy grid */
        if (_fission_rescaled) {
	    int lower_index = getEnergyGridIndex(energy);
	    double lower_xs = getFissionXS(lower_index);
	    double upper_xs = getFissionXS(lower_index + 1);
	    double delta_xs = upper_xs - lower_xs;
	    double slope = delta_xs / _delta_lethargy;
	    double lower_e = _start_lethargy + _delta_lethargy * lower_index;
	    double xs =  lower_xs + slope * (log10(energy) - lower_e);
	    log_printf(DEBUG, "low = %f, high = %f, xs = %f", lower_xs, 
			   upper_xs, xs);
	    return xs;
	}
	/* Use linear interpolation without a uniform lethargy grid */
	else {
	    return linearInterp<float, float, float>(_fission_xs_energies,
				      _fission_xs, _num_fission_xs, energy);
	}
    }
    /* If the isotope does not have a fission xs, return 0 */
    else
	return 0.0;
}


/**
 * @brief Returns the microscopic fission cross-section value for a certain 
 *        energy index in this isotope's rescaled uniform lethargy grid of 
 *        cross-section data.
 * @param energy_index the index into the energy array
 * @return the microscopic fission cross-section
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
 * @brief Returns the microscopic total cross-section value for some energy.
 * @details Uses linear interpolation to compute the cross-section at a
 *          a certain energy (eV).
 * @param energy the energy (eV) of interest
 * @return the microscopic total cross-section
 */
float Isotope::getTotalXS(float energy) const {

    /* If the total xs has been defined explicitly, use it to
     * linearly interpolate to find the total cross-section */
    if (_num_total_xs != 0) {
        /* Uses linear interpolation into a uniform lethargy grid */
        if (_rescaled) {
	    int lower_index = getEnergyGridIndex(energy);
	    double lower_xs = getTotalXS(lower_index);
	    double upper_xs = getTotalXS(lower_index + 1);
	    double delta_xs = upper_xs - lower_xs;
	    double slope = delta_xs / _delta_lethargy;
	    double lower_e = _start_lethargy + _delta_lethargy * lower_index;
	    double xs =  lower_xs + slope * (log10(energy) - lower_e);
	    log_printf(DEBUG, "low = %f, high = %f, xs = %f", lower_xs, 
			   upper_xs, xs);
	    return xs;
	}
	/* Uses linear interpolation without a uniform lethargy grid */
        else {
    	    return linearInterp<float, float, float>(_total_xs_energies,
				 _total_xs, _num_total_xs, energy);
	}
    }

    /* If a total cross-section array has not been computed for this isotope, 
     * loop over all xs which have been defined and add them to a total xs */
    else {
	float total_xs = 0;
	total_xs += getAbsorptionXS(energy);
	total_xs += getElasticXS(energy);
	return total_xs;
    }
}


/**
 * @brief Returns the microscopic total cross-section value for a certain 
 *        energy index in this isotope's rescaled uniform lethargy grid of 
 *        cross-section data.
 * @param energy_index the index into the energy array
 * @return the microscopic total cross-section
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
 * @brief Returns the microscopic transport cross-section value for some energy.
 * @details Uses linear interpolation to compute the cross-section at a
 *          a certain energy (eV). The transport cross-section corrects for
 *          anisotropic scattering and is computed as follows:
 *
 *          \f$ \sigma_{tr} = \sigma_t - \left<\mu\right>\sigma_s \f$
 *
 * @param energy the energy (eV) of interest
 * @return the microscopic transport cross-section
 */
float Isotope::getTransportXS(float energy) const {
    return (getTotalXS(energy) - _mu_avg * getElasticXS(energy));
}


/**
 * @brief Returns the microscopic transport cross-section value for a certain 
 *        energy index in this isotope's rescaled uniform lethargy grid of 
 *        cross-section data.
 * @details The transport cross-section corrects for anisotropic scattering
 *          using a linear approximation and is computed as follows:
 *
 *          \f$ \sigma_{tr} = \sigma_t - \left<\mu\right>\sigma_s \f$
 *
 * @param energy_index the index into the energy array
 * @return the microscopic transport cross-section
 */
float Isotope::getTransportXS(int energy_index) const {
    return (getTotalXS(energy_index) - _mu_avg * getElasticXS(energy_index));
}


/**
 * @brief This method returns true if the thermal scattering distributions
 *        for this isotope are to be used when sampling outgoing collision 
 *        energy.
 * @return boolean if the thermal scattering distributions exist
 */
bool Isotope::usesThermalScattering() {
    return _use_thermal_scattering;
}


/**
 * @brief This method returns whether or not the Isotope's cross-sections 
 *        have been rescaled to a uniform lethargy grid.
 * @return whether or not the cross-sections have been rescaled
 */
bool Isotope::isRescaled() const {
    return _rescaled;
}


/**
 * @brief Set the atomic number and update alpha, eta and rho.
 * @details Computes alpha, eta, rho and mu as follows:
 *
 *          \f$ \alpha = \left(\frac{A-1}{A+1}\right)^2 \f$
 *          \f$ \eta = \left(\frac{A+1}{2\sqrt{A}}\right)^2 \f$
 *          \f$ \rho = \left(\frac{A-1}{2\sqrt{A}}\right)^2 \f$
 *
 * @param A the isotope's atomic number
 */
void Isotope::setA(int A) {
    _A = A;
    _alpha = float(_A-1)/float(_A+1) * float(_A-1)/float(_A+1);
    _eta = (float(_A)+1.0) / (2.0 * sqrt(float(_A)));
    _rho = (float(_A)-1.0) / (2.0 * sqrt(float(_A)));
    _mu_avg = 2.0 / (3.0 * _A);
}



/**
 * @brief Sets the random number seed for reaction rate sampling.
 * @details This method is called by the Material::setRandomNumberSeed()
 *          method and should not be called directly by the user.
 * @param seed the random number seed
 */
void Isotope::setRandomNumberSeed(unsigned int seed) {
    _seed = seed;
}


/**
 * @brief Initializes the random number seed for random number sampling.
 * @details This method is called by the 
 *          Material::initializeRandomNumberGenerator()
 *          method and should not be called directly by the user.
 */
void Isotope::initializeRandomNumberGenerator() {
    srand(_seed);

    log_printf(NORMAL, "Initializing isotope %s random number seed to %d", _isotope_name, _seed);
}


/**
 * @brief Set the temperature of the isotope in degrees Kelvin.
 * @param T the temperature in degrees Kelvin
 */
void Isotope::setTemperature(float T) {
    _T = T;
}


/**
 * @brief Informs isotope not to use thermal scattering to sample outgoing
 *        collision energies.
 */
void Isotope::neglectThermalScattering() {
    _use_thermal_scattering = false;
}


/**
 * @brief Sets the thermal scattering high energy cutoff energy.
 * @param cutoff_energy the thermal scattering high energy cutoff energy (eV)
 */
void Isotope::setThermalScatteringCutoff(float cutoff_energy) {
    _thermal_cutoff = cutoff_energy;
}


/**
 * @brief Informs isotope to use thermal scattering to sample outgoing
 *        collision energies.
 */
void Isotope::useThermalScattering() {
    _use_thermal_scattering = true;
}


/**
 * @brief Inform isotope that it is fissionable.
 */
void Isotope::makeFissionable() {
    _fissionable = true;
}


/**
 * @brief Load the ENDF cross-section data from ASCII files into arrays
 *        for this isotope.
 * @details This method finds the appropriate ENDF data files for the isotope
 *          in the PINSPEC cross-section library based on the user-defined 
 *          name of the isotope. If the appropriate files are not found the
 *          method will return an exception. If only capture and elastic
 *          scattering cross-section data files are discovered in the 
 *          cross-section library then the isotope is not fissionable; otherwise
 *          if a fission cross-section file is found then the isotope is 
 *          fissionable. Finally, after all cross-sections are parsed in from
 *          data files, this method computes a total cross-section and an
 *          absorption cross-section and then rescales all cross-sections onto
 *          a uniform lethargy grid to allow for fast O(1) data lookup.
 */
void Isotope::loadXS() {

    log_printf(INFO, "Loading isotope %s", _isotope_name);

    /* initialize variables */
    std::string directory = getXSLibDirectory();
    std::string filename;
    struct stat buffer;
    float* energies;
    float* xs_values;

    /* Set this isotope's appropriate cross-sections using the data structure */
    /********************************** ELASTIC *******************************/
    /* Check whether an elastic cross-section file exists for isotope */

    filename = directory + _isotope_name + "-elastic.txt";

    if (stat(filename.c_str(), &buffer))
	log_printf(ERROR, "Unable to load elastic xs for isotope %s"
		   " since no data was found in the cross-section"
		   " file %s for this isotope", filename.c_str(), 
		   _isotope_name);

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


    /********************************** CAPTURE *******************************/
    /* Check whether a capture cross-section file exists for isotope */

    filename = directory + _isotope_name + "-capture.txt";

    if (stat(filename.c_str(), &buffer))
	log_printf(ERROR, "Unable to load capture xs for isotope %s"
		   " since no data was found in the cross-section"
		   " file %s for this isotope", filename.c_str(), _isotope_name); 

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

    /********************************* FISSION ********************************/
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
 * @brief Load the ENDF cross-section data for a particular cross-section from 
 *        an ASCII file into the appropriate array for this isotope.
 * @details This method finds the appropriate ENDF data file for the isotope
 *          in the PINSPEC cross-section library based on the user-defined name
 *          of the isotope as well as the type of cross-section input ('capture'
 *          'elastic', or 'fission'). If the appropriate file is not found the
 *          method will return an exception.  Finally, after the cross-section
 *          is parsed in from the data file, this method recomputes a total
 *          cross-section and an absorption cross-section and then rescales 
 *          all cross-sections onto a uniform lethargy grid to allow for 
 *          fast O(1) data lookup.
 * @param xs_type a character array for the cross-section type
 */
void Isotope::loadXS(char* xs_type) {

    std::string directory = getXSLibDirectory();
    std::string filename;
    struct stat buffer;
    float* energies;
    float* xs_values;

    /* Set this isotope's appropriate cross-section using the data structures */
    /********************************* ELASTIC ********************************/

    if (!strcmp(xs_type, "elastic")) {

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
    }


    /********************************** CAPTURE *******************************/
    if (!strcmp(xs_type, "capture")) {

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
    }


    /********************************* FISSION ********************************/
    if (!strcmp(xs_type, "fission")) {

        /* Check whether a fission cross-section file exists for isotope */
	filename = directory + _isotope_name + "-fission.txt";

	if (stat(filename.c_str(), &buffer))
	    log_printf(ERROR, "Unable to load fission xs for isotope %s"
		       " since no data was found in the cross-section"
		       " file %s for this isotope", 
                       filename.c_str(), _isotope_name); 

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
    
    /* Rescale the cross-sections onto a uniform lethargy grid */
    rescaleXS(pow(10., _start_lethargy), pow(10., _end_lethargy), _num_energies);

    return;
}



/**
 * @brief Set the elastic cross-section for this isotope.
 * @param elastic_xs a float array of microscopic elastic cross-sections
 * @param elastic_xs_energies a float array of energies (eV)
 * @param num_elastic_xs the number of elastic cross-section values
 */
void Isotope::setElasticXS(float* elastic_xs, float* elastic_xs_energies,
			   int num_elastic_xs) {
    _elastic_xs = elastic_xs;
    _elastic_xs_energies = elastic_xs_energies;
    _num_elastic_xs = num_elastic_xs;
    _elastic_rescaled = false;
}


/**
 * @brief Set the capture cross-section for this isotope.
 * @param capture_xs a float array of microscopic capture cross-sections
 * @param capture_xs_energies a float array of energies (eV)
 * @param num_capture_xs the number of capture cross-section
 */
void Isotope::setCaptureXS(float* capture_xs, float* capture_xs_energies,
			   int num_capture_xs) {
    _capture_xs = capture_xs;
    _capture_xs_energies = capture_xs_energies;
    _num_capture_xs = num_capture_xs;
    _capture_rescaled = false;
}


/**
 * @brief Set the fission cross-section for this isotope
 * @param fission_xs a float array of microscopic fission cross-sections
 * @param fission_xs_energies a float array of energies (eV)
 * @param num_fission_xs the number of fission cross-sections values
 */
void Isotope::setFissionXS(float* fission_xs, float* fission_xs_energies,
			   int num_fission_xs) {
    _fission_xs = fission_xs;
    _fission_xs_energies = fission_xs_energies;
    _num_fission_xs = num_fission_xs;
    _fission_rescaled = false;
}

/**
 * @brief Rescales all of the isotope's cross-sections onto a uniform
 *        lethargy grid.
 * @details Cross-section rescaling is useful because it allows for a fast
 *          O(1) table lookup (and linear interpolation) to compute 
 *          cross-section values for any given energy.
 * @param start_energy the highest lethargy value in the grid
 * @param end_energy the lowest lethargy value in the grid
 * @param num_energies the number of energies represented in the grid
 */
void Isotope::rescaleXS(float start_energy, float end_energy,
			int num_energies) {

    float* grid;
    float* new_energies;
    float* new_xs;

    _capture_rescaled = false;
    _elastic_rescaled = false;
    _fission_rescaled = false;

    /* Compute the uniform lethargy grid */
    grid = logspace<float, float>(start_energy, end_energy, num_energies);

    /* Capture xs */
    if (_num_capture_xs != 0) {

        new_energies = new float[num_energies];
	memcpy(new_energies, grid, sizeof(float)*num_energies);
	new_xs = new float[num_energies];

	for (int i=0; i < num_energies; i++)
	    new_xs[i] = getCaptureXS(new_energies[i]);

	_num_capture_xs = num_energies;
	delete [] _capture_xs_energies;
	delete [] _capture_xs;
	_capture_xs = new_xs;
	_capture_xs_energies = new_energies;
    }

    /* Elastic xs */
    if (_num_elastic_xs != 0) {
	new_energies = new float[num_energies];
	memcpy(new_energies, grid, sizeof(float)*num_energies);
	new_xs = new float[num_energies];

	for (int i=0; i < num_energies; i++)
	  new_xs[i] = getElasticXS(grid[i]);

	_num_elastic_xs = num_energies;
	delete [] _elastic_xs_energies;
	delete [] _elastic_xs;
	_elastic_xs = new_xs;
	_elastic_xs_energies = new_energies;
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
	delete [] _fission_xs;
	_fission_xs = new_xs;
	_fission_xs_energies = new_energies;
    }

    /* Compute absorption and total cross-sections */
    generateAbsorptionXS(start_energy, end_energy, num_energies);
    generateTotalXS(start_energy, end_energy, num_energies);

    /* Assign values for uniform lethargy grid parameters useful
     * for computing indices into the grid at a given energy */
    _start_lethargy = log10(start_energy);
    _end_lethargy = log10(end_energy);
    _delta_lethargy = (_end_lethargy - _start_lethargy) / _num_energies;
    _capture_rescaled = true;
    _elastic_rescaled = true;
    _fission_rescaled = true;

    delete [] grid;

    return;
}


/**
 * @brief Computes the microscopic absorption cross-section from 
 *        the isotope's capture and fission (if applicable) cross-sections.
 * @details This class method computes the absorption cross-section
 *          on a uniform lethargy grid.
 * @param start_energy the highest lethargy value in the grid
 * @param end_energy the lowest lethargy value in the grid
 * @param num_energies the number of energies represented in the grid
 */
void Isotope::generateAbsorptionXS(float start_energy, float end_energy, 
				   int num_energies) {

    /* Generate an absorption xs on the uniform energy/lethargy grid */
    float* new_energies = logspace<float, float>(start_energy, 
						 end_energy, num_energies);
    float* new_xs = new float[num_energies];

    _num_absorb_xs = num_energies;

    if (_fissionable)
        for (int i=0; i < num_energies; i++)
	    new_xs[i] = getCaptureXS(new_energies[i]) +                                                  getFissionXS(new_energies[i]);

    else
        for (int i=0; i < num_energies; i++)
	    new_xs[i] = getCaptureXS(new_energies[i]);


    _absorb_xs = new_xs;
    _absorb_xs_energies = new_energies;

    return;
}


/**
 * @brief Computes the microscopic total cross-section from the isotope's
 *        capture, elastic scatter and fission (if applicable) cross-sections.
 * @details This class method computes the total cross-section
 *          on a uniform lethargy grid.
 * @param start_energy the highest lethargy value in the grid
 * @param end_energy the lowest lethargy value in the grid
 * @param num_energies the number of energies represented in the grid
 */
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
 * @brief This method clones a given Isotope class object by executing a deep
 *        copy of all of the Isotope's class attributes and giving them to a new
 *        Isotope class object.
 * @return a pointer to the new cloned Isotope class object
 */
Isotope* Isotope::clone() {

    /* Allocate memory for the clone */
    Isotope* new_clone = new Isotope(_isotope_name);

    /* Set the clone's atomic number, number density */
    new_clone->setA(_A);
    new_clone->setTemperature(_T);

    /* Return a pointer to the cloned Isotope class */
    return new_clone;
}


/**
 * @brief Determines a random collision type based on the values of each of 
 *        the isotope's cross-section values at a given enery.
 * @param neutron a pointer to structure the of interest
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
 * @brief For a given neutron energy (eV) in a scattering collision, this
 *        function returns the outgoing energy in eV, \f$ E' \f$, for the 
 *        collision based on its thermal scattering distributions
 * @param energy the energy of the neutron of interest (eV)
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


/**
 * @brief This method initializes the probability distributions for thermal
 *        scattering.
 * @details It takes in arguments for the starting energy and end
 *          energy (ratios of kT) and the number of distributions which it
 *          uses to generate logarithmically spaced energies for the 
 *          distributions.
 * @param start_energy the first distribution's energy (ratio of kT)
 * @param end_energy the final distribution's energy (ratio of kT)
 * @param num_bins the number of bins per distribution
 * @param num_distributions the number of scattering distributions
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
			  end_energy/(_kB*_T), _num_thermal_cdfs);

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
 * @brief This function computes the thermal scattering probability for
 *        a ratio of initial to final energies.
 * @param E_prime_to_E a ratio of initial to final energies
 * @param dist_index the distribution of interest
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


/**
 * @brief Returns the number of thermal scattering CDFs.
 * @return the number of thermal scattering CDFs
 */
int Isotope::getNumThermalCDFs() {
    return _num_thermal_cdfs;
}


/**
 * @brief Returns the number of energy bins for each thermal scattering CDF.
 * @return the number of energy bins per thermal scattering CDF
 */
int Isotope::getNumThermalCDFBins() {
    return _num_thermal_cdf_bins;
}


/**
 * @brief Loads an input array with the values for each of the isotope's 
 *        thermal CDFs.
 * @details This method is intended to make the CDF data available to the
 *          PINSPEC user in Python. Although this function appears to require
 *          two input arguments - the cdfs array and the length of the array -
 *          in reality it only requires one argument for the array in Python.
 *          This method would be called in Python as follows:
 *
 * @code
 *          num_cdfs = isotope.getNumThermalCDFs()
 *          num_bins = isotope.getNumThermalCDFBins()
 *          cdfs = numpy.zeros(num_cdfs * num_bins)
 *          isotope.retrieveThermalCDFs(cdfs)
 * @endcode
 *
 * @param cdfs an input array for to fill with CDF values
 * @param num_values the number of CDF bins multiplied by the number of CDFs
 */
void Isotope::retrieveThermalCDFs(float* cdfs, int num_values) {
  
    for (int i=0; i < _num_thermal_cdfs; i++) {
        for (int j=0; j < _num_thermal_cdf_bins; j++)
            cdfs[i*_num_thermal_cdf_bins + j] = _thermal_cdfs[i][j];
    }  
}


/**
 * @brief Loads an input array with the energies for each of the isotope's 
 *        thermal PDFs.
 * @details This method is intended to make the PDF data available to the
 *          PINSPEC user in Python. Although this function appears to require
 *          two input arguments - the PDFs array and the length of the 
 *          array - in reality it only requires one argument for the array in 
 *          Python. This method would be called in Python as follows:
 *
 * @code
 *          num_cdfs = isotope.getNumThermalCDFs()
 *          num_bins = isotope.getNumThermalCDFBins()
 *          pdfs = numpy.zeros(num_cdfs * num_bins)
 *          isotope.retrieveThermalDistributions(pdfs)
 * @endcode
 *
 * @param pdfs an input array for to fill with CDF values
 * @param num_values the number of CDF bins multiplied by the number of CDFs
 */
void Isotope::retrieveThermalPDFs(float* pdfs, int num_values) {

    for (int i=0; i < _num_thermal_cdfs; i++) {
        for (int j=0; j < _num_thermal_cdf_bins; j++)
            pdfs[i*_num_thermal_cdf_bins + j] = 
	                        _thermal_dist[i*_num_thermal_cdf_bins + j];
    }  
}


/**
 * @brief Loads an input array with the \f$ \frac{E}{kT} \f$ values for each
 *        thermal scattering CDF. 
 * @details This method is intended to make data available to the PINSPEC
 *          user in Python. Although this function appears to require two
 *          input argument, in reality it only requires one argument for the
 *          array in Python. This method would be called in Python as follows:
 *
 * @code
 *          num_cdfs = isotope.getNumThermalCDFs()
 *          E_to_kT = numpy.zeros(num_cdfs)
 *          isotope.retrieveEtokT(E_to_kT)
 * @endcode
 *
 * @param E_to_kT an array of \f$ \frac{E}{kT} \f$ values
 * @param num_cdfs the number of thermal scattering CDFs
 */
void Isotope::retrieveEtokT(float* E_to_kT, int num_cdfs) {
    for (int i=0; i < _num_thermal_cdfs; i++)
        E_to_kT[i] = _E_to_kT[i];
}


/**
 * @brief Loads an input array with the \f$ \frac{E'}{E} \f$ values for each
 *        thermal scattering CDF.
 * @details This method is intended to make data available to the PINSPEC
 *          user in Python. Although this function appears to require two
 *          input arguments, in reality it only requires one argument for 
 *          the array in Python. This method would be called in Python as 
 *          follows:
 * @code
 *          num_bins = isotope.getNumThermalCDFBins()
 *          Eprime_to_E = numpy.zeros(num_bins)
 *          isotope.retrieveEprimeToE(E_prime_toE)
 * @endcode
 *
 * @param Eprime_to_E an array of \f$ \frac{E'}{E} \f$ values
 * @param num_bins the number of bins per thermal scattering CDF
 */
void Isotope::retrieveEprimeToE(float* Eprime_to_E, int num_bins) {

    for (int i=0; i < _num_thermal_cdf_bins; i++)
        Eprime_to_E[i] = _Eprime_to_E[i];
}



/**
 * @brief For a given neutron, this method samples a distance traveled to
 *        next collision.
 * @param neutron the neutron struct of interest
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
 * @brief For a given energy, this method samples a collision type, updates 
 *        the neutron's outgoing energy and kills the neutron if it sampled 
 *        absorption or leakage.
 * @param neutron the neutron struct of interest
 */
void Isotope::collideNeutron(neutron* neutron) {

    neutron->_old_energy = neutron->_energy;

    /* Obtain collision type and kill neutron if absorbed*/
    sampleCollisionType(neutron);

    /* Update the neutron's energy and direction on every collision,
     * including absorptions, to avoid having an if-else statement to 
     * check if it was absorbed or not. This is an optimization and
     * shouldn't slow the code down too much since most collisions will
     * result in scattering rather than capture or fission. */

    /* If the neutron is in an INFINITE_HOMOGENEOUS or HOMOGENEOUS_EQUIVALENCE 
     * geometry then it's azimuthal angle phi will have been initialized to 
     * -1 and we can simply uniformly sample an outgoing energy */
    if (neutron->_surface == NULL) {
        /* Sample outgoing energy uniformally between [alpha*E, E] */
        float alpha = getAlpha();
        double random = (float)(rand()) / (float)(RAND_MAX);

        /* Asymptotic elastic scattering above 4 eV */
        if (neutron->_energy > 4.0 || !_use_thermal_scattering)
    	    neutron->_energy *= (alpha + (1.0 - alpha) * random);
        else
    	    neutron->_energy = getThermalScatteringEnergy(neutron->_energy);
    }

    /* If the neutron is in a HETEROGENEOUS geometry then we must update
     * the azimuthal and polar angles for the neutron using isotropic
     * scattering in CM */
    else {
        /* Isotropic (in 2*pi) scattering for the azimuthal angle */
        float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
        float cos_phi = cos(phi);
        float sin_phi = sin(phi);

        /* Update neutron's direction vector - assume that scattering is
         * isotropic in center of mass for the polar angle */
        float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
        float mu_l = (1.0 + _A*mu_cm)/(sqrt(_A_squared + 2.0 * _A * mu_cm+1.0));
        float sqrt_mu_l_squared = sqrt(1.0 - mu_l * mu_l);

        float u = neutron->_u;
        float v = neutron->_v;
        float w = neutron->_w;
        float sqrt_w_squared = sqrt(1.0 - w * w);  
  
        neutron->_u = mu_l * u + ((sqrt_mu_l_squared * (u * w * cos_phi - 
					        v * sin_phi)) / sqrt_w_squared);
        neutron->_v = mu_l * v + ((sqrt_mu_l_squared * (v * w * cos_phi + 
					        u * sin_phi)) / sqrt_w_squared);
        neutron->_w = mu_l * w + (sqrt_mu_l_squared * sqrt_w_squared * cos_phi);

        float direction_norm = norm3D<float>(neutron->_u, neutron->_v, neutron->_w);
        neutron->_u /= direction_norm;
	neutron->_v /= direction_norm;
	neutron->_w /= direction_norm;

	//        neutron->_mu = cos(neutron->_w / 
	//		   norm2D<float>(neutron->_u, neutron->_v));
        //neutron->_phi = atan2(neutron->_v, neutron->_u);

        /* Correct for atan2's result in [-pi, pi] to interval [0, 2*pi] */
	//    if (neutron->_v <= 0.0)
	//  neutron->_phi += 2.0 * M_PI;
    
        /* Asymptotic elastic scattering above 4 eV */
        if (neutron->_energy > 4.0 || !_use_thermal_scattering)
            neutron->_energy *= (_A_squared + 2 *_A*mu_cm + 1.0) 
    	                    / _A_plus_one_squared;
        else
            neutron->_energy = getThermalScatteringEnergy(neutron->_energy);
    }

    neutron->_energy += 1E-7;  //FIXME: temp bug fix for zeroed out energy
    neutron->_collided = true;

    return;
}
