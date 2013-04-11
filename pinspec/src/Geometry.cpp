#include "Geometry.h"

/**
 * @brief Geomtery constructor.
 * @details Sets a default number of neutrons per batch (10,000), number of
 *          batches (10) and number of threads (1).
 */	
Geometry::Geometry(spatialType spatial_type) {

    /* Set defaults for the number of neutrons, batches, and threads */
    _num_neutrons_per_batch = 10000;
    _num_batches = 10;
    _num_threads = 1;

    _spatial_type = spatial_type;

    /* Initialize regions to null */
    _infinite_medium = NULL;
    _fuel = NULL;
    _moderator = NULL;

    /* Initialize a fissioner with a fission spectrum and sampling CDF */
    _fissioner = new Fissioner();
    
    /* Default for axial leakage is zero */
    _buckling_squared = 0.0;
}


/**
 * @brief Destructor deletes the Fissioner and lets SWIG delete 
 *        the regions, materials, isotopes and allies during 
 *        garbage collection.
 */
Geometry::~Geometry() {
    delete _fissioner;
}


/**
 * @brief Returns the number of neutrons per batch for this simulation.
 * @return the number of neutrons per batch
 */
int Geometry::getNumNeutronsPerBatch() {
    return _num_neutrons_per_batch;
}


/**
 * @brief Returns the total number of neutrons for this simulation.
 * @return the total number of neutrons for this simulation
 */
int Geometry::getTotalNumNeutrons() {
    return _num_neutrons_per_batch * _num_batches;
}


/**
 * @brief Returns the number of batches of neutrons for this simulation.
 * @return the number of batches
 */
int Geometry::getNumBatches() {
    return _num_batches;
}


/**
 * @brief Returns the number of parallel threads for this simulation.
 * @return the number of threads
 */
int Geometry::getNumThreads() {
    return _num_threads;
}


/**
 * @brief Return the spatial type of Geometry 
 *        (INFINITE_HOMOGENEOUS, HOMOGENEOUS_EQUIVALENCE or HETEROGENEOUS).
 * @return the spatial type
 */
spatialType Geometry::getSpatialType() {
    return _spatial_type;
}


/**
 * @brief Returns the square of the geometric buckling for this geometry.
 * @return the square of the geometric buckling
 */
float Geometry::getBucklingSquared() {
    return _buckling_squared;
}


/**
 * @brief Returns the total volume occuppied by the geometry.
 * @return The total volume for the geometry
 */
float Geometry::getVolume() {
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        return _infinite_medium->getVolume();
    else if (_spatial_type == HOMOGENEOUS_EQUIVALENCE)
        return _fuel->getVolume() + _moderator->getVolume();
    else
        return 1.0;     //FIXME: Update this for heterogeneous case when implemented
}



/**
 * @brief Sets the square of the geometric buckling for the geometry.
 * @param buckling_squared the square of the geometric buckling
 */ 
void Geometry::setBucklingSquared(float buckling_squared) {
    _buckling_squared = buckling_squared;
}


/**
 * @brief Sets the number of neutrons per batch for this simulation.
 * @param num_neutrons_per_batch the number of neutrons per batch
 */
void Geometry::setNeutronsPerBatch(int num_neutrons_per_batch) {
    _num_neutrons_per_batch = num_neutrons_per_batch;
}

/**
 * @brief Sets the number of batches for this simulation.
 * @param num_batches the number of batches
 */
void Geometry::setNumBatches(int num_batches) {
    _num_batches = num_batches;
}


/**
 * @brief Sets the number of batches for this simulation.
 * @param num_threads the number of batches	
 */
void Geometry::setNumThreads(int num_threads) {
    _num_threads = num_threads;
}


/**
 * @brief Sets the dancoff factor and computes the escape cross-section, 
 *        beta, alpha1 and alpha2 parameters used for a two region 
 *        heterogeneous-homogeneous pin cell simulation.
 * @param dancoff the dancoff factor
 */
void Geometry::setDancoffFactor(float dancoff) {

    _dancoff = dancoff;

    /* Two region homogeneous equivalence parameters */
    float A = (1.0 - dancoff) / dancoff;
    _sigma_e = 1.0 / (2.0 * _fuel->getFuelRadius());
    _alpha1 = ((5.0*A + 6.0) - sqrt(A*A + 36.0*A + 36.0)) / (2.0*(A+1.0));
    _alpha2 = ((5.0*A + 6.0) + sqrt(A*A + 36.0*A + 36.0)) / (2.0*(A+1.0));
    _beta = (((4.0*A + 6.0) / (A + 1.0)) - _alpha1) / (_alpha2 - _alpha1);
}


/**
 * @brief Set the geometry's spatial type
 *        (INFINITE_HOMOGENEOUS, HOMOGENEOUS_EQUIVALENCE or HETEROGENEOUS).
 * @param spatial_type the spatial type
 */
void Geometry::setSpatialType(spatialType spatial_type) {
    if (spatial_type == INFINITE_HOMOGENEOUS && _fuel != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " INFINITE_HOMOGENEOUS since it contains a FUEL Region %s",
		   _fuel->getRegionName());
    else if (spatial_type == INFINITE_HOMOGENEOUS && _moderator != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " INFINITE_HOMOGENEOUS since it contains a MODERATOR"
		   " Region %s", _moderator->getRegionName());
    else if (spatial_type == HOMOGENEOUS_EQUIVALENCE && 
	     _infinite_medium != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " HOMOGENEOUS_EQUIVALENCE since it contains an "
		   "INIFINTE Region %s", _infinite_medium->getRegionName());
    else if (spatial_type == HETEROGENEOUS && _infinite_medium != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " HETEROGENEOUS since it contains an"
		   " INFINITE_HOMOGENEOUS Region %s", 
		   _infinite_medium->getRegionName());
    else
      _spatial_type = spatial_type;
}


/**
 * @brief Adds a new region to the geometry.
 * @details Checks to make sure that the region type 
 *          (INFINITE, FUEL, MODERATOR) does not conflict with other regions
 *          that have already been added to the geometry
 * @param region the region to add to the geometry
 */
void Geometry::addRegion(Region* region) {

    if (region->getRegionType() == INFINITE) {
        if (_fuel != NULL)
	    log_printf(ERROR, "Unable to add an INFINITE type region %s"
		       " to the geometry since it contains a"
		       " FUEL type region %s", region->getRegionName(), 
		       _fuel->getRegionName());
	else if (_moderator != NULL)
	    log_printf(ERROR, "Unable to add an INFINITE type region %s"
		       " to the geometry since it contains a"
		       " MODERATOR type region %s", region->getRegionName(), 
		       _moderator->getRegionName());
	else if (_infinite_medium == NULL)
	    _infinite_medium = region;
	else
	    log_printf(ERROR, "Unable to add a second INFINITE type region %s"
		       " to the geometry since it already contains"
		       " region %s", region->getRegionName(), 
		       _infinite_medium->getRegionName());
    }

    else if (region->getRegionType() == FUEL) {
        if (_infinite_medium != NULL)
	    log_printf(ERROR, "Unable to add a FUEL type region %s"
		       " to the geometry since it contains an"
		       " INFINITE_HOMOGENEOUS type region %s", 
		       region->getRegionName(), 
		       _infinite_medium->getRegionName());
	else if (_fuel == NULL)
	    _fuel = region;
	else
	    log_printf(ERROR, "Unable to add a second FUEL type region %s"
		       " to the geometry since it already contains"
		       " region %s", region->getRegionName(), 
		       _fuel->getRegionName());
    }

    else if (region->getRegionType() == MODERATOR) {
        if (_infinite_medium != NULL)
	    log_printf(ERROR, "Unable to add a MODERATOR type region %s"
		       " to the geometry since it contains an"
		       " INFINITE_HOMOGENEOUS type region %s", 
		       region->getRegionName(), 
		       _infinite_medium->getRegionName());
	else if (_moderator == NULL)
	    _moderator = region;
	else
	    log_printf(ERROR, "Unable to add a second MODERATOR type region %s"
		       " to the geometry since it already contains"
		       " region %s", region->getRegionName(), 
		       _fuel->getRegionName());
    }
    else
        log_printf(ERROR, "Unable to add Region %s since it does not have a"
		   " Region type", region->getRegionName());
}


/**
 * @brief The primary Monte Carlo kernel for a PINSPEC simulation.
 * @details This method executes an appropriate Monte Carlo kernel depending
 *          on the geometry's spatial type. This method loops over batches
 *          and neutrons and collides each neutron in the appropriate region
 *          until it is absorbed, while tallying all user-specific quanties
 *          throughout.
 */
void Geometry::runMonteCarloSimulation() {

    if (_infinite_medium == NULL && _fuel == NULL && _moderator == NULL)
        log_printf(ERROR, "Unable to run Monte Carlo simulation since the"
		   " Geometry does not contain any Regions");

    int start_batch = 0;
    int end_batch = _num_batches;
    neutron* curr;
    bool precision_triggered = true;
    Timer timer;
    timer.start();
    TallyBank* tally_bank = TallyBank::Get();

    /* Print report to the screen */
        log_printf(TITLE, "Beginning PINSPEC Monte Carlo Simulation...");
    log_printf(NORMAL, "# neutrons / batch = %d     # batches = %d     "
                     "# threads = %d", _num_neutrons_per_batch, 
                        _num_batches, _num_threads);
    log_printf(SEPARATOR, "");

    omp_set_num_threads(_num_threads);
    tally_bank->initializeBatchTallies(_num_batches);

    /*************************************************************************/
    /************************   INFINITE_HOMOGENEOUS *************************/
    /*************************************************************************/

    /* If we are running an infinite medium spectral calculation */
    if (_spatial_type == INFINITE_HOMOGENEOUS){
        _infinite_medium->setBucklingSquared(_buckling_squared);

	while (precision_triggered) {

            #pragma omp parallel
            {
                #pragma omp for private(curr)
	        for (int i=start_batch; i < end_batch; i++) {

		    log_printf(INFO, "Thread %d/%d running batch %d", 
			       omp_get_thread_num()+1, 
			       omp_get_num_threads(), i);

		    curr = initializeNewNeutron();
	            curr->_batch_num = i;

		    for (int j=0; j < _num_neutrons_per_batch; j++) {

		        /* Initialize neutron energy [ev] from Watt spectrum */
		        curr->_energy = _fissioner->emitNeutroneV();
		        curr->_alive = true;
		        curr->_region = _infinite_medium;
			
			/* While the neutron is still alive, collide it. All
                         * tallying and collision physics take place within
			 * the Region, Material, and Isotope classes filling
                         * the Geometry */
			while (curr->_alive == true) {
			    _infinite_medium->collideNeutron(curr);
			    tally_bank->tally(curr);
			}
		    }
		}
            }

            tally_bank->computeScaledBatchStatistics(_num_neutrons_per_batch);
            if (!tally_bank->isPrecisionTriggered())
                precision_triggered = false;
            else {
                tally_bank->incrementNumBatches(_num_batches);
                start_batch = end_batch;
                end_batch += _num_batches;                
            }
        }
    }



    /*************************************************************************/
    /**********************   HOMOGENEOUS_EQUIVALENCE ************************/
    /*************************************************************************/

    /* If we are running homogeneous equivalence spectral calculation */
    else if (_spatial_type == HOMOGENEOUS_EQUIVALENCE) {
	
        /* Check that all necessary parameters have been set */
        if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
	    log_printf(ERROR, "Unable to run a HOMOGENEOUS_EQUIVALENCE type "
		       " simulation since beta, sigma_e, alpha1, or alpha2 for "
		       " have not yet been set for the geometry");

        else if (_fuel->getFuelRadius() == 0)
	    log_printf(ERROR, "Unable to run simulation since region %s "
		       " does not know the fuel radius", 
		       _fuel->getRegionName());
        else if (_moderator->getFuelRadius() == 0)
	    log_printf(ERROR, "Unable to run simulation since region %s "
		       " does not know the fuel radius", 
		       _moderator->getRegionName());
        else if (_fuel->getPitch() == 0)
	    log_printf(ERROR, "Unable to run simulation since region %s "
		       " does not know the pin pitch", _fuel->getRegionName());
        else if (_fuel->getPitch() == 0)
	    log_printf(ERROR, "Unable to run simulation since region %s "
		       " does not know the pin pitch", 
		       _moderator->getRegionName());

	_fuel->setBucklingSquared(_buckling_squared);
        _moderator->setBucklingSquared(_buckling_squared);
        initializeProbModFuelRatios();

	/* Initialize neutrons from fission spectrum for each thread */
	/* Loop over batches */
	/* Loop over neutrons per batch*/		
	float p_ff;
	float p_mf;
	float test;
       
        while (precision_triggered) {
            #pragma omp parallel
	    {
                #pragma omp for private(curr, p_ff, p_mf, test)
	        for (int i=start_batch; i < end_batch; i++) {
	  
	            log_printf(INFO, "Thread %d/%d running batch %d", 
	                      omp_get_thread_num()+1, omp_get_num_threads(), i);

		    curr = initializeNewNeutron();
		    curr->_batch_num = i;
				
		    for (int j=0; j < _num_neutrons_per_batch; j++) {

	            /* Initialize neutron's energy [ev] from Watt spectrum */
	            curr->_energy = _fissioner->emitNeutroneV();
		    curr->_alive = true;
		    curr->_region = _fuel;
			
		    /* While the neutron is still alive, collide it. All
		     * tallying and collision physics take place within
		     * the Region, Material, and Isotope classes filling
		     * the Geometry */
		    while (curr->_alive == true) {
	  
	                /* Determine if neutron collided in fuel or moderator */
	                p_ff = computeFuelFuelCollisionProb(curr);
	                p_mf = computeModeratorFuelCollisionProb(curr);
	                test = float(rand()) / RAND_MAX;
	  
	                /* If the neutron is in the fuel */
                        if (curr->_region == _fuel) {
	  
	                    /* If test is larger than p_ff, move to moderator */
	                    if (test > p_ff) {
	                        curr->_region = _moderator;
		                _moderator->collideNeutron(curr);
	                    }
	                    else
		                _fuel->collideNeutron(curr);
	                }
		        /* If the neutron is in the moderator */
			else {
	  
	                   /* If test is larger than p_mf, move to fuel */
	                   if (test < p_mf) {
	                       curr->_region = _fuel;
	                       _fuel->collideNeutron(curr);
	                   }
			   else
			       _moderator->collideNeutron(curr);
	                }

		        /* Tally the neutron collision */
			tally_bank->tally(curr);
	            }
	        }
	    }
	}
        tally_bank->computeScaledBatchStatistics(_num_neutrons_per_batch);
        if (!tally_bank->isPrecisionTriggered())
            precision_triggered = false;
        else {
	    tally_bank->incrementNumBatches(_num_batches);
	    start_batch = end_batch;
	    end_batch += _num_batches;                
	}
        }
    }



    /*************************************************************************/
    /***************************   HETEROGENEOUS *****************************/
    /*************************************************************************/

    //TODO: Truly heterogeneous geometry is NOT yet implemented!
    //The following code is simply to stub out what will become
    //the heterogeneous monte carlo kernel

    /* If we are running homogeneous equivalence spectral calculation */
    else if (_spatial_type == HETEROGENEOUS) {

    /* Check that all necessary parameters have been set */
    if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
        log_printf(ERROR, "Unable to run a HOMOGENEOUS_EQUIVALENCE type "
	    " simulation since beta, sigma_e, alpha1, or alpha2 for "
	    " have not yet been set for the geometry");

    else if (_fuel->getFuelRadius() == 0)
        log_printf(ERROR, "Unable to run simulation since region %s "
	    " does not know the fuel radius", _fuel->getRegionName());
    else if (_moderator->getFuelRadius() == 0)
        log_printf(ERROR, "Unable to run simulation since region %s "
	    " does not know the fuel radius", _moderator->getRegionName());
    else if (_fuel->getPitch() == 0)
        log_printf(ERROR, "Unable to run simulation since region %s "
	    " does not know the pin pitch", _fuel->getRegionName());
    else if (_fuel->getPitch() == 0)
        log_printf(ERROR, "Unable to run simulation since region %s "
	    " does not know the pin pitch", _moderator->getRegionName());

    /* Initialize neutrons from fission spectrum for each thread */
    /* Loop over batches */
    /* Loop over neutrons per batch*/		
    neutron* curr = initializeNewNeutron();

    while (precision_triggered) {
        for (int i=start_batch; i < end_batch; i++) {
	  
	    log_printf(INFO, "Thread %d/%d running batch %d", 
	    omp_get_thread_num()+1, omp_get_num_threads(), i);

	    curr->_batch_num = i;

	    for (int j=0; j < _num_neutrons_per_batch; j++) {

	    /* Initialize this neutron's energy [ev] from Watt spectrum */
	    curr->_energy = _fissioner->emitNeutroneV();
	    curr->_alive = true;
			
	    /* While the neutron is still alive, collide it. All
	     * tallying and collision physics take place within
	     * the Region, Material, and Isotope classes filling
	     * the Geometry */
	    while (curr->_alive == true) { 
	        /* Find the region the neutron is currently within */
	        /* Collide the neutron in the fuel or moderator */
	            _fuel->collideNeutron(curr);
	            tally_bank->tally(curr);
	    }
	}
	}
	tally_bank->computeScaledBatchStatistics(_num_neutrons_per_batch);
	if (!tally_bank->isPrecisionTriggered())
	    precision_triggered = false;
	else {
	  tally_bank->incrementNumBatches(_num_batches);
	  start_batch = end_batch;
	  end_batch += _num_batches;                
	}
        }


        /* Clone each Tally for each fuel ring */
        /* Clone each Tally for each moderator ring */
        /* Initialize neutrons from fission spectrum for each thread */
        /* Loop over batches */
        /* Loop over neutrons per batch*/
        /* Sample distance traveled in 3D */
        /* Compute new x, y coordinates */
        /* If neutron is in fuel, check if crossed the boundary of the fuel */
        /* If neutron is in moderator, check if it crossed into the fuel */
        /* If neutron stays current region, call region collideNeutron method */
	}


    /* Compute batch statistics for all Tallies in this simulation */
    timer.stop();
    log_printf(NORMAL, "PINSPEC simulated %.0f neutrons / sec in %f sec", 
        _num_neutrons_per_batch * end_batch / timer.getTime(),
	timer.getTime());
}


/**
 * @brief Initializes a pre-computed array of moderator to fuel first flight
 *        probability ratios.
 * @details The pre-computaiton of the ratios is an optimization to save time
 *          for the homogeneous-heterogeneous equivalence geometry type.
 */
void Geometry::initializeProbModFuelRatios() {

    Material* mod = _moderator->getMaterial();
    Material* fuel = _fuel->getMaterial();
    float v_mod = _moderator->getVolume();
    float v_fuel = _fuel->getVolume();
    _num_ratios = mod->getNumXSEnergies((char*)"elastic");

    /* Allocate memory for ratios */    
    _pmf_ratios = new float[_num_ratios];

    /* Set energy bounds and delta to allow for O(1) lookup of ratio */
    float* xs_energies = new float[_num_ratios];
    mod->retrieveXSEnergies(xs_energies, _num_ratios, (char*)"elastic");
    _start_energy = xs_energies[0];
    _end_energy = xs_energies[_num_ratios-1];

    if (_scale_type == EQUAL)
        _delta_energy = (xs_energies[_num_ratios-1] - 
	                xs_energies[0]) / _num_ratios;
    else {
        _start_energy = log10(xs_energies[0]);
	_end_energy = log10(xs_energies[_num_ratios-1]);
	_delta_energy = (_end_energy - _start_energy) / _num_ratios;
    }

    delete xs_energies;

    /* Loop over all xs energies and compute the P_mf ratios */
    for (int i=0; i < _num_ratios; i++)
        _pmf_ratios[i] = (fuel->getTotalMacroXS(i) * v_fuel) / 
                         (mod->getTotalMacroXS(i) * v_mod);

    return;
}


/**
 * @brief This function computes the two-region fuel-to-fuel collision 
 *        probability for a two-region pin cell simulation. It uses 
 *        Carlvik's two-term rational model.
 * @param neutron the neutron of interest
 * @return the fuel-to-fuel collision probability at the neutron's energy
 */
float Geometry::computeFuelFuelCollisionProb(neutron* neutron) {
    float energy = neutron->_energy;
    float p_ff;
    float sigma_tot_fuel = _fuel->getMaterial()->getTotalMacroXS(energy);
    p_ff = ((_beta*sigma_tot_fuel) / (_alpha1*_sigma_e + sigma_tot_fuel)) +
          ((1.0 - _beta)*sigma_tot_fuel / (_alpha2*_sigma_e + sigma_tot_fuel));
    log_printf(DEBUG, "sigma_tot_fuel = %f, p_ff = %f", sigma_tot_fuel, p_ff);
    return p_ff;
}


/**
 * @brief This function computes the two-region moderator-to-fuel collision
 *        probability for a two-region pin cell simulation. It uses Carlvik's
 *        two-term rational model.
 * @param neutron the neutron of interest
 * @return the moderator-to-fuel collision probability at the neutron's energy
 */
float Geometry::computeModeratorFuelCollisionProb(neutron* neutron) {
    float energy = neutron->_energy;
    float p_ff = computeFuelFuelCollisionProb(neutron);
    int index = getEnergyGridIndex(energy);
    float pmf_ratio = _pmf_ratios[index];
    float p_mf = (1.0 - p_ff) * pmf_ratio;
    return p_mf;
}
