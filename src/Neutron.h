/*
 * Neutron.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

//#include "Region.h"
class Region;

/* Structure to represent a neutron */
struct neutron {
	int _batch_num;
	float _x, _y, _z;
    float _u, _v, _w;
	float _mu, _phi;
	float _energy;
	bool _alive;
    Region* _region;
	bool _in_fuel;
};

neutron* initializeNewNeutron();


#endif /* NEUTRON_H_ */
