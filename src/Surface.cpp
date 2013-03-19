/*
 * Surface.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Surface.h"

/**
 * Surface constructor sets default values
 */
Surface::Surface() {
	_boundary_type = VACUUM;		/* Default boundary type */
}


/**
 * Surface destructor
 */
Surface::~Surface() { }



/**
 * Returns this Surface's boundary type (REFLECTIVE, VACUUM, or INTERFACE)
 * @return Surface boundary type
 */
boundaryType Surface::getBoundaryType() const {
    return _boundary_type;
}


/**
 * Sets the boundary type for this Surface (REFLECTIVE, VACUUM, or INTERFACE)
 * @param type the boundary type
 */
void Surface::setBoundaryType(boundaryType type) {
    _boundary_type = type;
}


/******************************************************************************
 *****************************   XPlane   *************************************
 *****************************************************************************/

XPlane::XPlane() {
	_x = 0.0;
};

XPlane::~XPlane() { };

float XPlane::getX() {
	return _x;
}

void XPlane::setX(float x) {
	_x = x;
}


float XPlane::computeNearestDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
    float x = neutron->_x;
    float u = neutron->_u;

	if ((x < _x && u > 0.0) || (x > _x && u < 0.0))
		dist = (_x - x) / u;

	return dist;
}


/**
 * Checks whether a certain x values is on the XPlane
 * (within numerical error)
 * @param x the value to check
 * @return if on the XPlane (true), otherwise false
 */
bool XPlane::onSurface(neutron* neutron) {
	if (fabs(_x - neutron->_x) < 1E-6)
		return true;

	return false;
}



/******************************************************************************
 *****************************   YPlane   *************************************
 *****************************************************************************/

YPlane::YPlane() {
	_y = 0.0;
};

YPlane::~YPlane() { };


float YPlane::getY() {
	return _y;
}

void YPlane::setY(float y) {
	_y = y;
}


float YPlane::computeNearestDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
    float y = neutron->_y;
    float v = neutron->_v;

	if ((y < _y && v > 0.0) || (y > _y && v < 0.0))
		dist = (_y - y) / v;

	return dist;
}


/**
 * Checks whether a certain y values is on the YPlane
 * (within numerical error)
 * @param x the value to check
 * @return if on the YPlane (true), otherwise false
 */
bool YPlane::onSurface(neutron* neutron) {
	if (fabs(_y - neutron->_y) < 1E-6)
		return true;

	return false;
}



/******************************************************************************
 *****************************   Circle   *************************************
 *****************************************************************************/

Circle::Circle() {
	_x0 = 0;
	_y0 = 0;
	_r = 0;
	_r_squared = 0;
};


Circle::~Circle() { };


float Circle::getX0() {
	return _x0;
}


float Circle::getY0() {
	return _y0;
}


float Circle::getRadius() {
	return _r;
}


void Circle::setX0(float x0) {
	_x0 = x0;
}


void Circle::setY0(float y0) {
	_y0 = y0;
}


void Circle::setRadius(float r) {
	_r = r;
	_r_squared = r*r;
}


float Circle::computeNearestDistance(neutron* neutron) {

	float x = neutron->_x;
	float y = neutron->_y;
	float u = neutron->_u;
	float v = neutron->_v;

	float a = u*u + v*v;
	float b = 2.0*x*u - 2.0*_x0*u + 2.0*y*v - 2.0*_y0*v;
	float c = x*x + _x0*_x0 - 2.0*_x0*x + y*y + _y0*_y0 - 2.0*_y0*y - _r_squared;

	float discr = b*b - 4.0*a*c;

	/* There is not an intersection point */
	if (discr < 0.0)
		return std::numeric_limits<float>::infinity();

	/* There is one intersection point */
	else if (discr == 0.0) {

		float dist = -b / (2.0*a);

		if (dist >= 0.0)
           return dist;
        else
           return std::numeric_limits<float>::infinity();
	}

	/* There are two intersection points */
	else {
        discr = sqrt(discr);
		float dist1 = (-b + discr) / (2.0*a);
		float dist2 = (-b - discr) / (2.0*a);

        if (dist1 <= dist1 && dist1 > 0.0)
            return dist1;
        else if (dist2 <= dist1 && dist2 > 0.0)
            return dist2;
        else
            return std::numeric_limits<float>::infinity();
    }
}


bool Circle::onSurface(neutron* neutron) {

	float r_squared = (neutron->_y - _y0) * (neutron->_y - _y0) + 
                        (neutron->_x - _x0) * (neutron->_x - _x0);

	if (fabs(_r_squared - r_squared) < 1E-6)
		return true;
	else
		return false;
}
