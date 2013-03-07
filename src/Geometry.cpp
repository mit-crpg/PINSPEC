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
 * Sets the escape cross-section, beta, alpha1 and alpha2 parameters
 * used for a two region pin cell simulation. The Geometry must have
 * HOMOGENEOUS_EQUIVALENCE spatial type set for this to work.
 * @param sigma_e the escape cross-section
 * @param beta Carlvik's beta parameter
 * @param alpha1 Carlvik's alpha1 parameter
 * @param alpha2 Carlvik's alpha2 parameter
 */
void Geometry::setTwoRegionPinCellParams(float sigma_e, float beta, 
								float alpha1, float alpha2) {
	if (_spatial_type == INFINITE_HOMOGENEOUS) {
		log_printf(ERROR, "Cannot set two region HOMOGENEOUS_EQUIVALENCE"
						" parameters for an INFINITE Geometry");
	}
	if (_spatial_type == HETEROGENEOUS) {
		log_printf(ERROR, "Cannot set two region HOMOGENEOUS_EQUIVALENCE"
						" parameters for a HETEROGENEOUS Geometry");
	}

	_sigma_e = sigma_e;
	_beta = beta;
	_alpha1 = alpha1;
	_alpha2 = alpha2;
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


void Geometry::runMonteCarloSimulation() {

	if (_infinite_medium == NULL && _fuel == NULL && _moderator == NULL)
		log_printf(ERROR, "Unable to run Monte Carlo simulation since the"
						" Geometry does not contain any Regions");

	/* Print report to the screen */
	log_printf(INFO, "Beginning PINSPEC Monte Carlo Simulation...");
	log_printf(INFO, "# neutrons / batch = %d\t\t# batches = %d\t\t# threads = %d",
						_num_neutrons_per_batch, _num_batches, _num_threads);


	/*************************************************************************/
	/************************   INFINITE_HOMOGENEOUS *************************/
	/*************************************************************************/

    /* If we are running an infinite medium spectral calculation */
	if (_spatial_type == INFINITE_HOMOGENEOUS){

        /* Inform Region of the number of batches to run - this is needed to
         * initialize the Tally classes for this many batches */
        _infinite_medium->setNumBatches(_num_batches);

		neutron* curr = initializeNewNeutron();

		for (int i=0; i < _num_batches; i++) {

			log_printf(INFO, "Batch #: %d", i);

			curr->_batch_num = i;

			for (int j=0; j < _num_neutrons_per_batch; j++) {

				/* Initialize this neutron's energy [ev] from Watt spectrum */
				curr->_energy = _fissioner->emitNeutroneV();
				curr->_alive = true;
			
				/* While the neutron is still alive, collide it. All
                 * tallying and collision physics take place within
				 * the Region, Material, and Isotope classes filling
                 * the Geometry
                 */
				while (curr->_alive == true) {
					_infinite_medium->collideNeutron(curr);
                    log_printf(DEBUG, "Updated neutron energy to %f", curr->_energy);
				}
			}
		}
	}



	/*************************************************************************/
	/**********************   HOMOGENEOUS_EQUIVALENCE ************************/
	/*************************************************************************/

	/* If we are running homogeneous equivalence spectral calculation */
	else if (_spatial_type == HOMOGENEOUS_EQUIVALENCE) {

        /* Inform Regions of the number of batches to run - this is needed to
         * initialize the Tally classes for this many batches */
        _fuel->setNumBatches(_num_batches);
        _moderator->setNumBatches(_num_batches);

		/* Check that all necessary parameters have been set */
		if (_beta <= 0 || _sigma_e <= 0 || _alpha1 <= 0 || _alpha2 <= 0)
			log_printf(ERROR, "Unable to run a HOMOGENEOUS_EQUIVALENCE type "
					" simulation since beta, sigma_e, alpha1, or alpha2 for "
					" have not yet been set for the geometry");

		/* Initialize neutrons from fission spectrum for each thread */
		/* Loop over batches */
		/* Loop over neutrons per batch*/		
		neutron* curr = initializeNewNeutron();
		float p_ff;
		float p_mf;
		float test;

		for (int i=0; i < _num_batches; i++) {

			log_printf(INFO, "Batch #: %d", i);

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
						if (test < p_mf)
							curr->_in_fuel = true;
					}

					/* Collide the neutron in the fuel or moderator */
					if (curr->_in_fuel)
						_fuel->collideNeutron(curr);
					else
						_moderator->collideNeutron(curr);
				}
			}
		}
    }



	/*************************************************************************/
	/***************************   HETEROGENEOUS *****************************/
	/*************************************************************************/

	/* If we are running homogeneous equivalence spectral calculation */
	else if (_spatial_type == HETEROGENEOUS) {

        /* Inform Regions of the number of batches to run - this is needed to
         * initialize the Tally classes for this many batches */
        _fuel->setNumBatches(_num_batches);
        _moderator->setNumBatches(_num_batches);

		/* Initialize neutrons from fission spectrum for each thread */
		/* Loop over batches */
		/* Loop over neutrons per batch*/
	}
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

    /* Call on the Region(s) inside the geometry to output statistics */
    if (_spatial_type == INFINITE_HOMOGENEOUS)
        _infinite_medium->outputBatchStatistics(directory, suffix);
    else {
        _fuel->outputBatchStatistics(directory, suffix);
        _moderator->outputBatchStatistics(directory, suffix);
    }

    return;   
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
	return p_ff;
}


/**
 * This function computes the two-region moderator-to-fuel collision
 * probability for a two-region pin cell simulation. It uses Carlvik's
 * two-term rational model.
 * @param energy the energy for a neutron in eV
 * @return the moderator-to-fuel collision probability at that energy
 */
float Geometry::computeModeratorFuelCollisionProb(float energy) {
	float p_mf;
	float p_ff = computeFuelFuelCollisionProb(energy);
	float p_fm = 1.0 - p_ff;
	float tot_sigma_f = _fuel->getMaterial()->getTotalMacroXS(energy);
	float tot_sigma_mod = _moderator->getMaterial()->getTotalMacroXS(energy);
	float v_mod = _moderator->getVolume();
	p_mf = p_fm*(tot_sigma_f*_fuel->getVolume()) / (tot_sigma_mod*v_mod);
	return p_mf;
}
