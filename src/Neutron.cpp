/*
 * Neutron.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Neutron.h"

/**
 * Creates a new neutron with a set of default attribute values
 * @return a pointer to the new neutron
 */
neutron* initializeNewNeutron(int batch_num) {

	/* Allocate memory for new neutron struct */
	neutron* new_neutron = new neutron;

	new_neutron->_batch_num = batch_num;

	/* Assign default values to neutron attributes */
	new_neutron->_x = 0.0;
	new_neutron->_y = 0.0;
	new_neutron->_z = 0.0;
	new_neutron->_energy = 0.0;
	new_neutron->_phi = 0.0;
	new_neutron->_mu = 0.0;
	new_neutron->_weight = 1.0;
	new_neutron->_alive = true;

	return new_neutron;
}

