#include "Neutron.h"

/**
 * @brief Creates a new neutron with a set of default attribute values.
 * @details Intializes the neutron with 0.0 for its position and direction
 *          vector, energy, and batch number and makes the neutron's
 *          alive attribute true.
 * @return a pointer to the new neutron
 */
neutron* initializeNewNeutron() {

    /* Allocate memory for new neutron struct */
    neutron* new_neutron = new neutron;

    /* Assign default values to neutron attributes */
    new_neutron->_x = 0.0;
    new_neutron->_y = 0.0;
    new_neutron->_z = 0.0;

    new_neutron->_u = 0.0;
    new_neutron->_v = 0.0;
    new_neutron->_w = 0.0;

    new_neutron->_phi = 0.0;
    new_neutron->_mu = 0.0;
    new_neutron->_energy = 0.0;
    new_neutron->_alive = true;
    new_neutron->_batch_num = 0.0;

    new_neutron->_old_energy = 0.0;
    new_neutron->_total_xs = 0.0;

    return new_neutron;
}
