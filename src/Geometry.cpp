/*
 * Region.cpp
 *
 *  Created on: Mar 6, 2013
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Geometry.h"


/**
 * Geomtery constructor sets empty default material name
 */	
Geometry::Geometry() {

	/* Set defaults for the number of neutrons, batches, and threads */
	_num_neutrons_per_batch = 10000;
	_num_batches = 10;
	_num_threads = 1;

	_spatial_type = INFINITE_HOMOGENEOUS;

	/* Initialize regions to null */
	_infinite_medium = NULL;
	_fuel = NULL;
	_moderator = NULL;

	/* Initialize a fissioner with a fission spectrum and sampling CDF */
	_fissioner = new Fissioner();
	_fissioner->setNumBins(10000);
	_fissioner->setEMax(20);
	_fissioner->buildCDF();

}


/**
 * Geometry destructor deletes the Fissioner and lets SWIG delete 
 * the Regions, Materials, Isotopes and Tallies
 */
Geometry::~Geometry() {

	delete _fissioner;
}


/**
 * Returns the number of neutrons per batch for this simulation
 * @return the number of neutrons per batch
 */
int Geometry::getNumNeutronsPerBatch() {
	return _num_neutrons_per_batch;
}


/**
 * Returns the total number of neutrons for this simulation
 * @return the total number of neutrons for this simulation
 */
int Geometry::getTotalNumNeutrons() {
	return _num_neutrons_per_batch * _num_batches;
}


/**
 * Returns the number of batches of neutrons for this simulation
 * @return the number of batches
 */
int Geometry::getNumBatches() {
	return _num_batches;
}


/**
 * Returns the number of parallel threads for this simulation
 * @return the number of threads
 */
int Geometry::getNumThreads() {
	return _num_threads;
}


/**
 * Return the spatial type of Geometry 
 * (INFINITE_HOMOGENEOUS, HOMOGENEOUS_EQUIVALENCE or HETEROGENEOUS)
 * @return the spatial type
 */
spatialType Geometry::getSpatialType() {
	return _spatial_type;
}


float Geometry::getTotalMacroXS(float energy, Region* region) {
    return region->getTotalMacroXS(energy);
}


float Geometry::getTotalMacroXS(int energy_index, Region* region) {
    return region->getTotalMacroXS(energy_index);
}


float Geometry::getTotalMicroXS(float energy, Region* region) {
    return region->getTotalMicroXS(energy);
}


float Geometry::getTotalMicroXS(int energy_index, Region* region) {
    return region->getTotalMicroXS(energy_index);
}


float Geometry::getElasticMacroXS(float energy, Region* region) {
    return region->getTotalMacroXS(energy);
}


float Geometry::getElasticMacroXS(int energy_index, Region* region) {
    return region->getTotalMacroXS(energy_index);
}


float Geometry::getElasticMicroXS(float energy, Region* region) {
    return region->getTotalMicroXS(energy);
}


float Geometry::getElasticMicroXS(int energy_index, Region* region) {
    return region->getTotalMicroXS(energy_index);
}


float Geometry::getAbsorptionMacroXS(float energy, Region* region) {
    return region->getTotalMacroXS(energy);
}


float Geometry::getAbsorptionMacroXS(int energy_index, Region* region) {
    return region->getTotalMacroXS(energy_index);
}


float Geometry::getAbsorptionMicroXS(float energy, Region* region) {
    return region->getTotalMicroXS(energy);
}


float Geometry::getAbsorptionMicroXS(int energy_index, Region* region) {
    return region->getTotalMicroXS(energy_index);
}


float Geometry::getCaptureMacroXS(float energy, Region* region) {
    return region->getTotalMacroXS(energy);
}


float Geometry::getCaptureMacroXS(int energy_index, Region* region) {
    return region->getTotalMacroXS(energy_index);
}


float Geometry::getCaptureMicroXS(float energy, Region* region) {
    return region->getTotalMicroXS(energy);
}


float Geometry::getCaptureMicroXS(int energy_index, Region* region) {
    return region->getTotalMicroXS(energy_index);
}


float Geometry::getFissionMacroXS(float energy, Region* region) {
    return region->getTotalMacroXS(energy);
}


float Geometry::getFissionMacroXS(int energy_index, Region* region) {
    return region->getTotalMacroXS(energy_index);
}


float Geometry::getFissionMicroXS(float energy, Region* region) {
    return region->getTotalMicroXS(energy);
}


float Geometry::getFissionMicroXS(int energy_index, Region* region) {
    return region->getTotalMicroXS(energy_index);
}


float Geometry::getTransportMacroXS(float energy, Region* region) {
    return region->getTotalMacroXS(energy);
}


float Geometry::getTransportMacroXS(int energy_index, Region* region) {
    return region->getTotalMacroXS(energy_index);
}


float Geometry::getTransportMicroXS(float energy, Region* region) {
    return region->getTotalMicroXS(energy);
}


float Geometry::getTransportMicroXS(int energy_index, Region* region) {
    return region->getTotalMicroXS(energy_index);
}


/**
 * Sets the number of neutrons per batch for this simulation
 * @param the number of neutrons per batch
 */
void Geometry::setNeutronsPerBatch(int num_neutrons_per_batch) {
	_num_neutrons_per_batch = num_neutrons_per_batch;
}



/**
 * Sets the number of batches for this simulation
 * @param the number of batches
 */
void Geometry::setNumBatches(int num_batches) {
	_num_batches = num_batches;
}


/**
 * Sets the number of batches for this simulation
 * @param the number of batches	
 */
void Geometry::setNumThreads(int num_threads) {
	_num_threads = num_threads;
}


/**
 * Sets the dancoff factor and computes the escape cross-section, 
 * beta, alpha1 and alpha2 parameters used for a two region pin 
 * cell simulation.
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
 * Set the Geometry's spatial type 
 * (INFINITE_HOMOGENEOUS, HOMOGENEOUS_EQUIVALENCE or HETEROGENEOUS)
 * @param spatial_type the spatial type
 */
void Geometry::setSpatialType(spatialType spatial_type) {
	if (spatial_type == INFINITE_HOMOGENEOUS && _fuel != NULL)
		log_printf(ERROR, "Cannot set the geometry spatial type to be"
						" INFINITE_HOMOGENEOUS since it contains a FUEL Region %s",
						_fuel->getRegionName());
	else if (spatial_type == INFINITE_HOMOGENEOUS && _moderator != NULL)
		log_printf(ERROR, "Cannot set the geometry spatial type to be"
						" INFINITE_HOMOGENEOUS since it contains a MODERATOR Region %s",
						_moderator->getRegionName());
	else if (spatial_type == HOMOGENEOUS_EQUIVALENCE && _infinite_medium != NULL)
		log_printf(ERROR, "Cannot set the geometry spatial type to be"
						" HOMOGENEOUS_EQUIVALENCE since it contains an "
						"INIFINTE Region %s", 
						_infinite_medium->getRegionName());
	else if (spatial_type == HETEROGENEOUS && _infinite_medium != NULL)
		log_printf(ERROR, "Cannot set the geometry spatial type to be"
						" HETEROGENEOUS since it contains an"
						" INFINITE_HOMOGENEOUS Region %s", 
						_infinite_medium->getRegionName());
	else
		_spatial_type = spatial_type;
}


/**
 * Adds a new Region to the Geometry. Checks to make sure that the Region
 * type (INFINITE, FUEL, MODERATOR) does not conflict with other Regions
 * that have already been added to the Geometry
 * @param a region to add to the Geometry
 */
void Geometry::addRegion(Region* region) {

	if (region->getRegionType() == INFINITE) {
		if (_fuel != NULL)
			log_printf(ERROR, "Unable to add an INFINITE type region %s"
								" to the geometry since it contains a"
								" FUEL type region %s",
								region->getRegionName(), 
								_fuel->getRegionName());
		else if (_moderator != NULL)
			log_printf(ERROR, "Unable to add an INFINITE type region %s"
								" to the geometry since it contains a"
								" MODERATOR type region %s", 
								region->getRegionName(), 
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
 * Add a tally class object to  this Geometry's vector of Tallies
 * @param tally a pointer to a Tally object
 */
void Geometry::addTally(Tally *tally) {

    if (tally->getTallyDomainType() != GEOMETRY)
        log_printf(ERROR, "Unable to add Tally %s to Geometry since the Tally"
                        " is not for an GEOMETRY tally domain", 
                                        tally->getTallyName());
    _tallies.push_back(tally);
    return;
}


void Geometry::runMonteCarloSimulation() {

	if (_infinite_medium == NULL && _fuel == NULL && _moderator == NULL)
		log_printf(ERROR, "Unable to run Monte Carlo simulation since the"
						" Geometry does not contain any Regions");

    int start_batch = 0;
    int end_batch = _num_batches;
    neutron* curr;
    float sample;
    collisionType type;
    bool precision_triggered = true;
    Timer timer;
    timer.start();

	/* Print report to the screen */
	log_printf(NORMAL, "Beginning PINSPEC Monte Carlo Simulation...");
	log_printf(NORMAL, "# neutrons / batch = %d     # batches = %d     "
                      "# threads = %d", _num_neutrons_per_batch, 
                        _num_batches, _num_threads);

    
	omp_set_num_threads(_num_threads);
    initializeBatchTallies();

	/*************************************************************************/
	/************************   INFINITE_HOMOGENEOUS *************************/
	/*************************************************************************/


    /* If we are running an infinite medium spectral calculation */
    if (_spatial_type == INFINITE_HOMOGENEOUS){

        while (precision_triggered) {

            #pragma omp parallel
            {
	            #pragma omp for private(curr, sample, type)
	        for (int i=start_batch; i < end_batch; i++) {

		        log_printf(NORMAL, "Thread %d/%d running batch %d", 
                                        omp_get_thread_num()+1, 
                                        omp_get_num_threads(), i);

		        curr = initializeNewNeutron();
	            curr->_batch_num = i;

		        for (int j=0; j < _num_neutrons_per_batch; j++) {

			        /* Initialize neutron energy [ev] from Watt spectrum */
			        curr->_energy = _fissioner->emitNeutroneV();
			        curr->_alive = true;
			
			        /* While the neutron is still alive, collide it. All
                     * tallying and collision physics take place within
		   	         * the Region, Material, and Isotope classes filling
                     * the Geometry
                     */
			        while (curr->_alive == true) {
                        sample = curr->_energy;
				        type = _infinite_medium->collideNeutron(curr);
                        tally(sample, curr->_batch_num, 
                                            _infinite_medium, type);
			        }
		        }
            }
            }

            computeScaledBatchStatistics();
            if (!isPrecisionTriggered())
                precision_triggered = false;
            else {
                incrementNumBatches(_num_batches);
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
		float p_ff;
		float p_mf;
		float test;

        initializePmfRatios();

        while (precision_triggered) {
            #pragma omp parallel
            {
               #pragma omp for private(curr, p_ff, p_mf, test, sample, type)
		        for (int i=start_batch; i < end_batch; i++) {

		            log_printf(NORMAL, "Thread %d/%d running batch %d", 
                                  omp_get_thread_num()+1, omp_get_num_threads(), i);

	                curr = initializeNewNeutron();
			        curr->_batch_num = i;

			        for (int j=0; j < _num_neutrons_per_batch; j++) {

				        /* Initialize neutron's energy [ev] from Watt spectrum */
				        curr->_energy = _fissioner->emitNeutroneV();
				        curr->_alive = true;
				        curr->_in_fuel = true;
			
				        /* While the neutron is still alive, collide it. All
                         * tallying and collision physics take place within
				         * the Region, Material, and Isotope classes filling
                         * the Geometry
                         */
				        while (curr->_alive == true) {

					        /* Determine if neutron collided in fuel or moderator */
					        p_ff = computeFuelFuelCollisionProb(curr->_energy);
					        p_mf = computeModeratorFuelCollisionProb(curr->_energy);
					        test = float(rand()) / RAND_MAX;

					        /* If the neutron is in the fuel */
					        if (curr->_in_fuel) {

						        /* If test is larger than p_ff, move to moderator */
						        if (test > p_ff)
							        curr->_in_fuel = false;
					        }

					        /* If the neutron is in the moderator */
					        else {

						        /* If test is larger than p_mf, move to fuel */
						        if (test > p_mf)
							        curr->_in_fuel = true;
					        }

					        /* Collide the neutron in the fuel or moderator */
					        if (curr->_in_fuel) {
                                sample = curr->_energy;
						        type = _fuel->collideNeutron(curr);
                                tally(sample, curr->_batch_num, _fuel, type);
                            }
					        else {
                                sample = curr->_energy;
						        type = _moderator->collideNeutron(curr);
                                tally(sample, curr->_batch_num, _moderator, type);
                            }
				        }
			        }
		        }
            }
            computeScaledBatchStatistics();
            if (!isPrecisionTriggered())
                precision_triggered = false;
            else {
                incrementNumBatches(_num_batches);
                start_batch = end_batch;
                end_batch += _num_batches;                
            }
        }
    }



	/*************************************************************************/
	/***************************   HETEROGENEOUS *****************************/
	/*************************************************************************/

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
//		float p_ff;
//		float p_mf;
//		float test;

        while (precision_triggered) {
		    for (int i=start_batch; i < end_batch; i++) {

			        log_printf(NORMAL, "Thread %d/%d running batch %d", 
                                  omp_get_thread_num()+1, omp_get_num_threads(), i);

			    curr->_batch_num = i;

			    for (int j=0; j < _num_neutrons_per_batch; j++) {

				    /* Initialize this neutron's energy [ev] from Watt spectrum */
				    curr->_energy = _fissioner->emitNeutroneV();
				    curr->_alive = true;
				    curr->_in_fuel = true;
			
				    /* While the neutron is still alive, collide it. All
                     * tallying and collision physics take place within
				     * the Region, Material, and Isotope classes filling
                     * the Geometry
                     */
				    while (curr->_alive == true) { 
				        /* Collide the neutron in the fuel or moderator */
				        if (curr->_in_fuel) {
                            sample = curr->_energy;
					        type = _fuel->collideNeutron(curr);
                            tally(sample, curr->_batch_num, _fuel, type);
                        }
				        else {
                            sample = curr->_energy;
					        type = _moderator->collideNeutron(curr);
                            tally(sample, curr->_batch_num, _moderator, type);
                        }
                    }
                }
            }
            computeScaledBatchStatistics();
            if (!isPrecisionTriggered())
                precision_triggered = false;
            else {
                incrementNumBatches(_num_batches);
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
            /* If neutron is in fuel, check if it crossed the boundary of the fuel */
            /* If neutron is in moderator, check if it crossed into the fuel, or one of the cell boundaries */
            /* If neutron stays in current region, call region collideNeutron method */
            /* Need to update Region collideNeutron method to do additional tallying for each ring */
	}


    /* Compute batch statistics for all Tallies in this simulation */
    timer.stop();
    log_printf(NORMAL, "PINSPEC simulated %.0f neutrons / sec in %f sec", 
                _num_neutrons_per_batch * end_batch / timer.getTime(),
                timer.getTime());
}


bool Geometry::isPrecisionTriggered() {

    /* Check if the Geometry has any Tallies with a trigger on the precision */
    std::vector<Tally*>::iterator iter;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter) {
        if ((*iter)->isPrecisionTriggered())
            return true;
    }

    /* Check if the Regions have any Tallies with a trigger on the precision */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        return _infinite_medium->isPrecisionTriggered();
    else {
        return (_fuel->isPrecisionTriggered() 
                || _moderator->isPrecisionTriggered());
    }
}


void Geometry::computeBatchStatistics() {

    /* Compute statistics for each of this Geometry's Tallies */
    std::vector<Tally*>::iterator iter;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter)
        (*iter)->computeBatchStatistics();

    /* Compute statistics for Tallies inside the Region's inside the Geometry */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        _infinite_medium->computeBatchStatistics();
    else {
        _fuel->computeBatchStatistics();
        _moderator->computeBatchStatistics();
    }

    return;
}


void Geometry::computeScaledBatchStatistics() {

    /* Compute statistics for each of this Geometry's Tallies */
    std::vector<Tally*>::iterator iter;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter)
        (*iter)->computeScaledBatchStatistics(_num_neutrons_per_batch);

    /* Compute statistics for Tallies inside the Region's inside the Geometry */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        _infinite_medium->computeScaledBatchStatistics(_num_neutrons_per_batch);
    else {
        _fuel->computeScaledBatchStatistics(_num_neutrons_per_batch);
        _moderator->computeScaledBatchStatistics(_num_neutrons_per_batch);
    }

    return;
}


/**
 * Calls each of the Tally class objects in the simulation to output
 * their tallies and statistics to output files. If a user asks to
 * output the files to a directory which does not exist, this method
 * will create the directory.
 * @param directory the directory to write batch statistics files
 * @param suffix a string to attach to the end of each filename
 */
void Geometry::outputBatchStatistics(char* directory,  char* suffix) {

    /* Check to see if directory exists - if not, create it */
    struct stat st;
    if (!stat(directory, &st) == 0) {
		mkdir(directory, S_IRWXU);
    }

    /* Compute statistics for each of this Geometry's Tallies */
    std::vector<Tally*>::iterator iter;
    std::string filename;

	for (iter = _tallies.begin(); iter != _tallies.end(); ++iter) {
        filename = std::string(directory) + "/" + (*iter)->getTallyName()
                                               + "_" + suffix + ".txt";
        (*iter)->outputBatchStatistics(filename.c_str());
    }

    /* Compute statistics for Tallies inside the Region's inside the Geometry */

    /* Call on the Region(s) inside the geometry to output statistics */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        _infinite_medium->outputBatchStatistics(directory, suffix);
    else {
        _fuel->outputBatchStatistics(directory, suffix);
        _moderator->outputBatchStatistics(directory, suffix);
    }

    return;   
}


void Geometry::tally(float sample, int batch_num, 
                                        Region* region, collisionType type) {

	float total_xs = getTotalMacroXS(sample, region);
    std::vector<Tally*>::iterator iter;

    /* Tallies the event into the appropriate tally classes  */
    for (iter = _tallies.begin(); iter != _tallies.end(); iter ++) {
        Tally *tally = *iter;
        tallyType tally_type = tally->getTallyType();
        switch (tally_type) {
        case FLUX:
		    tally->weightedTally(sample, 1.0 / total_xs, batch_num);
        case COLLISION_RATE:
	        tally->weightedTally(sample, 1.0, batch_num);
        case ELASTIC_RATE:
	    if (type == ELASTIC)
	        tally->weightedTally(sample, 
		    getElasticMacroXS(sample, region) / total_xs, batch_num);
        case ABSORPTION_RATE:
//	    if (type == CAPTURE || type == FISSION)
	        tally->weightedTally(sample, 
		    getAbsorptionMacroXS(sample, region) / total_xs, batch_num);
        case CAPTURE_RATE:
	    if (type == CAPTURE)
	        tally->weightedTally(sample, 
            getCaptureMacroXS(sample, region) / total_xs, batch_num);
        case FISSION_RATE:
	    if (type == FISSION)
	        tally->weightedTally(sample, 
		    getFissionMacroXS(sample, region) / total_xs, batch_num);
        case TRANSPORT_RATE:
	    if (type == ELASTIC)
	        tally->weightedTally(sample, 
		    getTransportMacroXS(sample, region) / total_xs, batch_num);
        case DIFFUSION_RATE: /* FIXME */
	    if (type == ELASTIC)
	        tally->weightedTally(sample, 
		     1.0 / (3.0 * getTransportMacroXS(sample, region) * total_xs), batch_num); 
        case LEAKAGE_RATE:; /* FIXME */
        }
    }
}


void Geometry::initializeBatchTallies() {

    std::vector<Tally*>::iterator iter;

    /* Set the number of batches for all of this Geometry's Tallies */
    for (iter = _tallies.begin(); iter != _tallies.end(); iter ++)
        (*iter)->setNumBatches(_num_batches);

    /* Inform Region of the number of batches to run - this is needed to
     * initialize the Tally classes for this many batches */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        _infinite_medium->setNumBatches(_num_batches);
    else {
        _fuel->setNumBatches(_num_batches);
        _moderator->setNumBatches(_num_batches);
    }

}


void Geometry::initializePmfRatios() {

    Material* mod = _moderator->getMaterial();
    Material* fuel = _fuel->getMaterial();
    float v_mod = _moderator->getVolume();
    float v_fuel = _fuel->getVolume();
    _num_ratios = mod->getNumXSEnergies();
    _scale_type = mod->getEnergyGridScaleType();

    if (_scale_type == LOGARITHMIC)
        log_printf(NORMAL, "Scale type set to LOGARITHMIC for %d pmf ratios", _num_ratios);

    /* Allocate memory for ratios */    
    _pmf_ratios = new float[_num_ratios];

    /* Set energy bounds and delta to allow for O(1) lookup of ratio */
    float* xs_energies = new float[_num_ratios];
    mod->retrieveXSEnergies(xs_energies, _num_ratios);
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


void Geometry::incrementNumBatches(int num_batches) {

    std::vector<Tally*>::iterator iter;

    /* Update the number of batches for all of this Geometry's Tallies */
    for (iter = _tallies.begin(); iter != _tallies.end(); iter ++)
        (*iter)->incrementNumBatches(_num_batches);

    /* Inform Region of the number of batches to run - this is needed to
     * initialize the Tally classes for this many batches */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        _infinite_medium->incrementNumBatches(_num_batches);
    else {
        _fuel->incrementNumBatches(_num_batches);
        _moderator->incrementNumBatches(_num_batches);
    }
    
}


/**
 * This function computes the two-region fuel-to-fuel collision probability for
 * a two-region pin cell simulation. It uses Carlvik's two-term rational model.
 * @param energy the energy for a neutron in eV
 * @return the fuel-to-fuel collision probability at that energy
 */
float Geometry::computeFuelFuelCollisionProb(float energy) {
	float p_ff;
	float sigma_tot_fuel = _fuel->getMaterial()->getTotalMacroXS(energy);
	p_ff = ((_beta*sigma_tot_fuel) / (_alpha1*_sigma_e + sigma_tot_fuel)) +
		((1.0 - _beta)*sigma_tot_fuel / (_alpha2*_sigma_e + sigma_tot_fuel));
    log_printf(DEBUG, "sigma_tot_fuel = %f, p_ff = %f", sigma_tot_fuel, p_ff);
	return p_ff;
}


/**
 * This function computes the two-region moderator-to-fuel collision
 * probability for a two-region pin cell simulation. It uses Carlvik's
 * two-term rational model.
 * @param ene rgy the energy for a neutron in eV
 * @return the moderator-to-fuel collision probability at that energy
 */
float Geometry::computeModeratorFuelCollisionProb(float energy) {
	float p_ff = computeFuelFuelCollisionProb(energy);
    int index = getEnergyGridIndex(energy);
    float pmf_ratio = _pmf_ratios[index];       // FIXME: Use liner interpolation for more accuracy
    float p_mf = (1.0 - p_ff) * pmf_ratio;
    log_printf(DEBUG, "p_ff = %f, p_mf = %f", p_ff, p_mf);
	return p_mf;
}
