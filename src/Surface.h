/*
 * Surface.h
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <limits>
#include <math.h>
#include "Neutron.h"

#define PI_OVER_TWO 1.57079633
#define THREE_PI_OVER_TWO 4.71238898
#define TWO_PI 6.28318531
#define TINY_MOVE 1E-3


/* The types of boundaries */
typedef enum boundaryTypes {
    REFLECTIVE,
	VACUUM,
	INTERFACE
} boundaryType;


class Surface
{
protected:
	boundaryType _boundary_type;
public:
	Surface();
	virtual ~Surface();

    boundaryType getBoundaryType() const;

    void setBoundaryType(boundaryType type);

    virtual float computeNearestDistance(neutron* neutron) =0;
    virtual bool onSurface(neutron* neutron) =0;
};


class XPlane: public Surface {
private:
	float _x;
public:
	XPlane();
	virtual ~XPlane();
	float getX();
	void setX(float x);
    float computeNearestDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
};


class YPlane: public Surface {
private:
	float _y;
public:
	YPlane();
	virtual ~YPlane();
	float getY();
	void setY(float y);
    float computeNearestDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
};


class Circle: public Surface {
protected:
	float _r, _r_squared;
	float _x0, _y0;
public:
	Circle();
	virtual ~Circle();
	float getX0();
	float getY0();
	float getRadius();
	void setX0(float x0);
	void setY0(float y0);
	void setRadius(float r);
    float computeNearestDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
};


#endif /* SURFACE_H_ */
