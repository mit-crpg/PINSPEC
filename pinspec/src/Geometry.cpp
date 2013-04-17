#include "Geometry.h"

/**
 * @brief Geomtery constructor.
 * @details Sets a default number of neutrons per batch (10,000), number of
 *          batches (10) and number of threads (1). Sets the source sampling
 *          radius to 10 cm by default.
 */	
Geometry::Geometry(spatialType spatial_type, const char* name) {

    _geometry_name = name;

    /* Set defaults for the number of neutrons, batches, and threads */
    _num_neutrons_per_batch = 10000;
    _num_batches = 10;
    _num_threads = 1;

    _spatial_type = spatial_type;

    /* Initialize regions to null */
    _infinite_medium = NULL;
    _fuel = NULL;
    _moderator = NULL;

    /* Default equivalent geometry parameters */
    _fuel_radius = 0.45;
    _pitch = 1.26;

    /* The default number of first flight collision probabilities is zero */
    _num_prob = 0;

    /* Default for axial leakage is zero */
    _buckling_squared = 0.0;
   
    /* Default dancoff factor is non-physical to allow for error checking */
    _dancoff = -1.0;

    /* Initialize a fissioner with a fission spectrum and sampling CDF */
    _fissioner = new Fissioner();
    _source_sampling_radius = 2.0;
}


/**
 * @brief Destructor lets SWIG delete regions, materials, isotopes and
 *        tallies during garbage collection.
 */
Geometry::~Geometry() { 
    delete _fissioner;

    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE && _num_prob > 0) {
        delete [] _prob_ff;
        delete [] _prob_mf;
        delete [] _prob_energies;
    }
}


/**
 * @brief Returns the name of the geometry.
 * @return the name of the goemetry
 */
const char* Geometry::getName() {
    return _geometry_name;
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
 * @brief Returns the source sampling radius used for rejection sampling of
 *        random fission emission source sites.
 * @return the source sampling radius (cm)
 */
float Geometry::getSourceSamplingRadius() {
    return _source_sampling_radius;
}


/**
 * @brief Sets the name of the geometry.
 * @param name the name of the geometry
 */
void Geometry::setName(const char* name) {
    _geometry_name = name;
}


/**
 * @brief Sets the source sampling radius used for rejection sampling of
 *        random fission emission source sites.
 * @param radius the source sampling radius (cm)
 */
void Geometry::setSourceSamplingRadius(float radius) {
    _source_sampling_radius = radius;
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
 * @brief Sets the fuel pin radius.
 * @param radius the fuel pin radius (cm)
 */
void Geometry::setFuelPinRadius(float radius) {
    _fuel_radius = radius;
}


/**
 * @brief Sets the lattice pin cell pitch.
 * @param pitch the pin cell pitch (cm)
 */
void Geometry::setPinCellPitch(float pitch) {
  _pitch = pitch;
}


/**
 * @brief Sets the dancoff factor and computes the escape cross-section, 
 *        beta, alpha1 and alpha2 parameters used for a two region 
 *        heterogeneous-homogeneous pin cell simulation.
 * @param dancoff the dancoff factor
 */
void Geometry::setDancoffFactor(float dancoff) {
    if (dancoff > 1.0)
        log_printf(ERROR, "Unable to set a dancoff factor of %f since it is "
	       "greater than 1.0", dancoff);
    else if (dancoff < 0.0)
        log_printf(ERROR, "Unable to set a dancoff factor of %f since it is "
	       "less than 0.0", dancoff);
    else
        _dancoff = dancoff;
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
		   _fuel->getName());
    else if (spatial_type == INFINITE_HOMOGENEOUS && _moderator != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " INFINITE_HOMOGENEOUS since it contains a MODERATOR"
		   " Region %s", _moderator->getName());
    else if (spatial_type == HOMOGENEOUS_EQUIVALENCE && 
	     _infinite_medium != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " HOMOGENEOUS_EQUIVALENCE since it contains an "
		   "INIFINTE Region %s", _infinite_medium->getName());
    else if (spatial_type == HETEROGENEOUS && _infinite_medium != NULL)
        log_printf(ERROR, "Cannot set the geometry spatial type to be"
		   " HETEROGENEOUS since it contains an"
		   " INFINITE_HOMOGENEOUS Region %s", 
		   _infinite_medium->getName());
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

    if (region->getRegionType() == INFINITE_MEDIUM) {
        if (_fuel != NULL)
	    log_printf(ERROR, "Unable to add an INFINITE type region %s"
		       " to the geometry since it contains a"
		       " FUEL type region %s", region->getName(), 
		       _fuel->getName());
	else if (_moderator != NULL)
	    log_printf(ERROR, "Unable to add an INFINITE type region %s"
		       " to the geometry since it contains a"
		       " MODERATOR type region %s", region->getName(), 
		       _moderator->getName());
	else if (_infinite_medium == NULL)
	    _infinite_medium = static_cast<InfiniteMediumRegion*>(region);
	else
	    log_printf(ERROR, "Unable to add a second INFINITE type region %s"
		       " to the geometry since it already contains"
		       " region %s", region->getName(), 
		       _infinite_medium->getName());
    }

    else if (region->getRegionType() == EQUIVALENT_FUEL) {
        if (_infinite_medium != NULL)
	    log_printf(ERROR, "Unable to add a FUEL type region %s"
		       " to the geometry since it contains an"
		       " INFINITE_HOMOGENEOUS type region %s", 
		       region->getName(), 
		       _infinite_medium->getName());
	else if (_fuel == NULL)
	    _fuel = static_cast<EquivalenceRegion*>(region);
	else
	    log_printf(ERROR, "Unable to add a second FUEL type region %s"
		       " to the geometry since it already contains"
		       " region %s", region->getName(), 
		       _fuel->getName());
    }

    else if (region->getRegionType() == EQUIVALENT_MODERATOR) {
        if (_infinite_medium != NULL)
	    log_printf(ERROR, "Unable to add a MODERATOR type region %s"
		       " to the geometry since it contains an"
		       " INFINITE_HOMOGENEOUS type region %s", 
		       region->getName(), 
		       _infinite_medium->getName());
	else if (_moderator == NULL)
	  _moderator = static_cast<EquivalenceRegion*>(region);
	else
	    log_printf(ERROR, "Unable to add a second MODERATOR type region %s"
		       " to the geometry since it already contains"
		       " region %s", region->getName(), 
		       _fuel->getName());
    }
    else if (region->getRegionType() == BOUNDED_GENERAL) {
        _regions.push_back(static_cast<BoundedRegion*>(region));
    }
    else if (region->getRegionType() == BOUNDED_MODERATOR) {
        _regions.push_back(static_cast<BoundedRegion*>(region));
    }
    else if (region->getRegionType() == BOUNDED_FUEL) {
        _regions.push_back(static_cast<BoundedRegion*>(region));
    }
    else
        log_printf(ERROR, "Unable to add Region %s since it does not have a"
		   " Region type", region->getName());
}


/**
 * @brief Determines whether or not the geometry contains this neutron.
 * @details If the geometry contains this neutron's location, sets the 
 *          neutron struct's region pointer to this region and returns
 *          true. If no region is found, returns false.
 *
 * @param neutron the neutron of interest
 */
bool Geometry::contains(neutron* neutron) {

    /* Simpy return true if the geometry is infinite or uses 
     * heterogeneous-homogeneous equivalence theory */
    if (_spatial_type == INFINITE_HOMOGENEOUS 
        || _spatial_type == HOMOGENEOUS_EQUIVALENCE)
        return true;

    /* If the geometry is heterogeneous, loop over all regions to find
     * if any contains this neutron */
    std::vector<BoundedRegion*>::iterator iter;
    for (iter = _regions.begin(); iter != _regions.end(); ++iter) {
        if ((*iter)->contains(neutron)) {
	    neutron->_region = (*iter);
            return true;
        }
    }

    /* If no containing region was found, return false */
    return false;
}


/**
 * @brief Determine whether or not this point is contained in the geometry.
 * @details Returns true for INFINITE_HOMOGENEOUS and HOMOGENEOUS_EQUIVALENCE
 *          geometries. Returns true for HETEROGENEOUS geometries if the
 *          point is contained within one of the geomtry's regions.
 * @param x the x-coordinate of interest
 * @param y the y-coordinate of interest
 * @param z the z-coordinate of interest
 */
bool Geometry::contains(float x, float y, float z) {

    if (_spatial_type == INFINITE_HOMOGENEOUS)
        return true;
    else if (_spatial_type == HOMOGENEOUS_EQUIVALENCE)
      return true;
    else {

        /* Loop over all of a HETEROGENOUS geometery's regions */
        std::vector<BoundedRegion*>::iterator iter;
        for (iter = _regions.begin(); iter != _regions.end(); ++iter) {

	    /* If this region contains the neutron */
	    if ((*iter)->contains(x, y, z))
	      return true;
	}
        
        /* None of the geomtry's regions contained the neutron */
        return NULL;
    }
}


/**
 * @brief Finds the region containing a neutron.
 * @details Finds the region and sets the neutron struct's region pointer 
 *          to this region. If no region is found, throws exception.
 *
 * @param neutron the neutron of interest
 */
void Geometry::findContainingRegion(neutron* neutron) {

    /* Simpy return if the geometry is infinite or uses 
     * heterogeneous-homogeneous equivalence theory */
    if (_spatial_type == INFINITE_HOMOGENEOUS 
        || _spatial_type == HOMOGENEOUS_EQUIVALENCE)
            return;

    /* If the geometry is heterogeneous, loop over all regions to find
     * one which contains this neutron */
    std::vector<BoundedRegion*>::iterator iter;
    for (iter = _regions.begin(); iter != _regions.end(); ++iter) {
        if ((*iter)->contains(neutron)) {
	    neutron->_region = (*iter);
            return;
        }
    }

    /* If no containing region was found, throw exception */
    log_printf(ERROR, "Unable to find the region containing neutron at "
	       "(x,y,z) = (%f,%f,%f)", neutron->_x, neutron->_y, neutron->_z);
}

/**
 * @brief Finds the region containing a neutron.
 * @details Finds the region containing a neutron and returns a pointer to it
 *          or NULL if no region was found. Also returns NULL for 
 *          INFINITE_HOMOGENEOUS and HOMOGENEOUS_EQUIVALENCE geometry types.
 * @param x the x-coordinate of interest
 * @param y the y-coordinate of interest
 * @param z the z-coordinate of interest
 */
Region* Geometry::findContainingRegion(float x, float y, float z) {

    /* Return NULL for INFINITE_HOMOGENEOUS and HOMOGENEOUS_EQUIVALENCE */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        return NULL;
    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE)
        return NULL;

    /* If the geometry is HETEROGENEOUS, loop over all regions to find
     * one which contains this neutron */
    std::vector<BoundedRegion*>::iterator iter;
    for (iter = _regions.begin(); iter != _regions.end(); ++iter) {
        if ((*iter)->contains(x, y, z))
	    return (*iter);
    }

    return NULL;
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

    /*************************************************************************/
    /****************************  ERROR CHECKING  ***************************/
    /*************************************************************************/
    if (_spatial_type == INFINITE_HOMOGENEOUS && _infinite_medium == NULL)
        log_printf(ERROR, "Unable to run Monte Carlo simulation for an "
		   "INFINITE_HOMOGENEOUS geometry since it does not contain "
		   " an INFINITE_MEDIUM type region.");
    if (_spatial_type == INFINITE_HOMOGENEOUS && 
	  _infinite_medium->getMaterial() == NULL)
              log_printf(ERROR, "Unable to run Monte Carlo simulation for an "
		   "INFINITE_HOMOGENEOUS geometry since the infinite medium "
                   "region does not contain a material.");
    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE && _fuel == NULL)
        log_printf(ERROR, "Unable to run Monte Carlo simulation for a "
		   "HOMOGENEOUS_EQUIVALENCE geometry since it does not contain "
		   "an EQUIVALENT_FUEL type region.");
    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE && _moderator == NULL)
        log_printf(ERROR, "Unable to run Monte Carlo simulation for a "
		   "HOMOGENEOUS_EQUIVALENCE geometry since it does not contain "
		   "an EQUIVALENT_MODERATOR type region.");    
    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE && 
	  _fuel->getMaterial() == NULL) 
              log_printf(ERROR, "Unable to run Monte Carlo simulation for a "
		   "HOMOGENEOUS_EQUIVALENCE geometry since the fuel "
                   "region does not contain a material.");
    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE && 
	  _moderator->getMaterial() == NULL) 
              log_printf(ERROR, "Unable to run Monte Carlo simulation for a "
		   "HOMOGENEOUS_EQUIVALENCE geometry since the moderator "
                   "region does not contain a material.");
    if (_spatial_type == HOMOGENEOUS_EQUIVALENCE && _dancoff == -1.0)
            log_printf(ERROR, "Unable to run a HOMOGENEOUS_EQUIVALENCE type "
		       "simulation since the dancoff factor has not yet "
		       "been set for the geometry");
    if (_spatial_type == HETEROGENEOUS && _regions.size() == 0)
        log_printf(ERROR, "Unable to run a HETEROGENEOUS type simulation since "
                   " the geometry does not contain any BOUNDED type regions");


    /*************************************************************************/
    /****************************  SIMULATION SETUP  *************************/
    /*************************************************************************/
    int start_batch = 0;
    int end_batch = _num_batches;
    neutron curr;
    bool precision_triggered = true;
    Timer timer;
    timer.start();
    TallyBank* tally_bank = TallyBank::Get();

    if (_spatial_type == INFINITE_HOMOGENEOUS) {
        _infinite_medium->setBucklingSquared(_buckling_squared);
    }
    else if (_spatial_type == HOMOGENEOUS_EQUIVALENCE) {
        _fuel->setBucklingSquared(_buckling_squared);
	_moderator->setBucklingSquared(_buckling_squared);

        _fuel->setOtherRegion(_moderator);
        _moderator->setOtherRegion(_fuel);

	_fuel->setFuelPinRadius(_fuel_radius);
	_moderator->setFuelPinRadius(_fuel_radius);
	
	_fuel->setPinCellPitch(_pitch);
	_moderator->setPinCellPitch(_pitch);

        initializeProbModFuelRatios();
    }
    else {
        std::vector<BoundedRegion*>::iterator iter;
        for (iter = _regions.begin(); iter != _regions.end(); ++iter)
              (*iter)->setBucklingSquared(_buckling_squared);
    }

    tally_bank->initializeBatchTallies(_num_batches);

    omp_set_num_threads(_num_threads);

    /* Print report to the screen */
    log_printf(TITLE, "Beginning PINSPEC Monte Carlo Simulation...");
    log_printf(NORMAL, "# neutrons / batch = %d     # batches = %d     "
                     "# threads = %d", _num_neutrons_per_batch, 
                        _num_batches, _num_threads);
    log_printf(SEPARATOR, "");


    /*************************************************************************/
    /*************************   MONTE CARLO KERNEL   ************************/
    /*************************************************************************/

    while (precision_triggered) {

        #pragma omp parallel
        {
            #pragma omp for private(curr)
            for (int i=start_batch; i < end_batch; i++) {

	        log_printf(INFO, "Thread %d/%d running batch %d", 
                       omp_get_thread_num()+1, omp_get_num_threads(), i);

                curr._batch_num = i;

                for (int j=0; j < _num_neutrons_per_batch; j++) {

	            /* Initialize new source neutron */
                    initializeSourceNeutron(&curr);
                    
                    /* While the neutron is still alive, collide it. All
                     * tallying and collision physics take place within
		     * the region, material, and isotope classes filling
                     * the geometry */
                     while (curr._alive == true) {
                         findContainingRegion(&curr);
                         curr._region->collideNeutron(&curr);
                         tally_bank->tally(&curr);
                    }
                }
	    }
	}

        /* Compute batch statistics for all tallies in this simulation */
        tally_bank->computeScaledBatchStatistics(_num_neutrons_per_batch);
        if (!tally_bank->isPrecisionTriggered())
            precision_triggered = false;
        else {
            tally_bank->incrementNumBatches(_num_batches);
            start_batch = end_batch;
            end_batch += _num_batches;                
        }
    }

    timer.stop();

    log_printf(NORMAL, "PINSPEC simulated %.0f neutrons / sec in %f sec", 
	    _num_neutrons_per_batch * end_batch / timer.getTime(),
	    timer.getTime());

    return;
}	


/**
 * @brief Initializes a pre-computed array of moderator to first flight
 *        collsion probabilities using Carlvik's two term rational model.
 * @details The pre-computaiton of the probabilities is an optimization to 
 *          save time in the monte carlo kernel for the 
 *          homogeneous-heterogeneous equivalence geometry type. Each of the
 *          fuel and moderator equivalence theory regions contains a reference
 *          to the arrays of first flight collision probabilities in addition
 *          to the geometery. Everyone (both regions and the geometry)
 *          reference the same arrays to optimize cache performance. The
 *          geometry is in charge of deleting the memory for the arrays at
 *          the end of the simulation.
 */
void Geometry::initializeProbModFuelRatios() {

    /* Two region homogeneous equivalence parameters */
    float A = (1.0 - _dancoff) / _dancoff;
    _sigma_e = 1.0 / (2.0 * _fuel->getFuelPinRadius());
    _alpha1 = ((5.0*A + 6.0) - sqrt(A*A + 36.0*A + 36.0)) / (2.0*(A+1.0));
    _alpha2 = ((5.0*A + 6.0) + sqrt(A*A + 36.0*A + 36.0)) / (2.0*(A+1.0));
    _beta = (((4.0*A + 6.0) / (A + 1.0)) - _alpha1) / (_alpha2 - _alpha1);

    if (_spatial_type != HOMOGENEOUS_EQUIVALENCE)
        log_printf(ERROR, "Unable to initialize first flight collision "
            "probabilities for the geometry since it is not a "
	    "HOMOGENEOUS_EQUIVALENCE type geometry.");

    Material* mod = _moderator->getMaterial();
    Material* fuel = _fuel->getMaterial();
    float v_mod = _moderator->getVolume();
    float v_fuel = _fuel->getVolume();
    float sigma_tot_fuel;

    _num_prob = mod->getNumXSEnergies((char*)"elastic");

    /* Allocate memory for first flight collision probabilities */    
    float* prob_mf_ratios = new float[_num_prob];
    _prob_ff = new float[_num_prob];
    _prob_mf = new float[_num_prob];

    /* Set energy bounds and delta to allow for O(1) lookup of probabilities */
    _prob_energies = new float[_num_prob];
    mod->retrieveXSEnergies(_prob_energies, _num_prob, (char*)"elastic");

    /* Loop over all xs energies and compute the P_mf ratios */
    for (int i=0; i < _num_prob; i++) {
        prob_mf_ratios[i] = (fuel->getTotalMacroXS(i) * v_fuel) / 
                         (mod->getTotalMacroXS(i) * v_mod);

        sigma_tot_fuel = fuel->getTotalMacroXS(i);
        _prob_ff[i] = ((_beta * sigma_tot_fuel) / (_alpha1 *_sigma_e + 
		sigma_tot_fuel)) + ((1.0 - _beta) * sigma_tot_fuel / 
		(_alpha2 * _sigma_e + sigma_tot_fuel));
        _prob_mf[i] = (1.0 - _prob_ff[i]) * prob_mf_ratios[i];
    }

    /* Load first flight collision probabilities into each homogeneous
     * equivalent region */
    _moderator->setFirstFlightCollProb(_prob_ff, _prob_mf, 
		                       _prob_energies, _num_prob);
    _fuel->setFirstFlightCollProb(_prob_ff, _prob_mf, 
		                       _prob_energies, _num_prob);

    delete prob_mf_ratios;

    return;
}


/**
 * @brief Initializes a new source neutron within the geometry.
 * @details A source neutron initialized within the geometry will 
 *          have an energy (eV) from a Watt spectrum, an _alive attribute
 *          set to true, a _collided attribute set to false, and its
 *          _region pointer set to this region. The _material and _isotope 
 *          attributes will be set to NULL. For HETEROGENEOUS geometries, 
 *          this method uses rejection sampling to initialize the neutrons 
 *          location to a source site within the geometry with a non-zero 
 *          fission cross-section. An isotropic (in lab) direction vector is 
 *          sampled for the neutron's trajectory in 3D.
 * @param neutron the neutron of interest
 */
void Geometry::initializeSourceNeutron(neutron* neutron) {

    neutron->_energy = _fissioner->emitNeutroneV();
    neutron->_old_energy = neutron->_energy;
    neutron->_collided = false;
    neutron->_total_xs = 0.0;
    neutron->_path_length = 0.0;
    neutron->_alive = true;
    neutron->_material = NULL;
    neutron->_isotope = NULL;
    neutron->_surface = NULL;

    if (_spatial_type == INFINITE_HOMOGENEOUS) {
        neutron->_region = _infinite_medium;
        neutron->_surface = NULL;
    }

    else if (_spatial_type == HOMOGENEOUS_EQUIVALENCE) {
        neutron->_region = _fuel;
        neutron->_surface = NULL;
    }

    /* Use rejection sampling to sample a source site within a HETEROGENEOUS 
     * geometry with a non-zero fission cross-section */
    else {
        int i;

        for (i=0; i < 1000; i++) {

	    /* Uniformly sample within a sphere */
            /* Azimuthal angle in xy-plane */
	    float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;

            /* Polar angle with respect to z-axis */
            float cos_theta = (float(rand()) / RAND_MAX) * 2.0 - 1.0; 
            float theta = acos(cos_theta);

            /* Radius from sphere center */
            float u = float(rand()) / RAND_MAX; 
            float radius = _source_sampling_radius * pow(u, (1./3.));

            /* Spherical to cartesian coordinate conversion */
            neutron->_x = radius * sin(theta) * cos(phi);
            neutron->_y = radius * sin(theta) * sin(phi);
            neutron->_z = radius * cos(theta);
     
            /* If we successfully found a source site within the geometry with
  	     * a non-zero thermal fission cross-section, break the loop */
            if (contains(neutron) && 
	        neutron->_region->getFissionMacroXS(0.0253f) > 0.0)
	            break;
        }

        if (i == 1000)
            log_printf(ERROR, "Unable to sample a source site with rejection "
	         "sampling using a sampling radius of %f cm after 1000 "
	         "attempts.", _source_sampling_radius);

        /* Randomly sample a direction vector */
        neutron->_u = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
        neutron->_v = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
        neutron->_w = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
    }

    return;
}
