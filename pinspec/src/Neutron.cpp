#include "Neutron.h"

/**
 * @brief Create a new empty neutron struct.
 * @return the pointer to the new neutron struct
 */
neutron* createNewNeutron() {
    neutron* neut = new neutron;

    neut->_batch_num = 0;
    neut->_energy = 0.0;
    neut->_old_energy = 0.0;

    return neut;
}
