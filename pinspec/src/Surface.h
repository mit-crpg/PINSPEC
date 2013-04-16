/**
 * @file Surface.h
 * @brief The Surface class
 * @author William Boyd (wboyd@mit.edu)
 * @date March 15, 2012
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <limits>
#include <math.h>
#include "vector.h"
#include "Neutron.h"

/** \f$ \frac{\pi}{2}\f$ */
#define PI_OVER_TWO 1.57079633
/** \f$ \frac{3\pi}{2}\f$ */
#define THREE_PI_OVER_TWO 4.71238898
/** \f$ 2\pi \f$ */
#define TWO_PI 6.28318531
/** A small move to push a particle across a Surface into a new Region */
#define TINY_MOVE 1E-5


/**
 * @enum surfaceTypes
 * @brief Bounding surface types
 */

/**
 * @var surfaceType
 * @brief Bounding surface type
 */
typedef enum surfaceTypes {
    /** A plane perpendicular to the x-axis */
    XPLANE,
    /** A plane perpendicular to the y-axis */
    YPLANE,
    /** An infinite cylinder parallel to the z-axis */
    ZCYLINDER
} surfaceType;


/**
 * @enum boundaryTypes
 * @brief Boundary condition types
 */

/**
 * @var boundaryType
 * @brief Boundary condition type
 */
typedef enum boundaryTypes {
    /** Reflective boundary conditions */ 
    REFLECTIVE,
    /** Vacuum boundary conditions */
    VACUUM,
    /** Interface boundary conditions */
    INTERFACE
} boundaryType;


/**
 * @class Surface Surface.h "pinspec/src/Surface.h"
 * @brief The Surface represents a quadratic surface in the xy-plane.
 * @details The Surface is an abstract class with stub for the function
 *          calls needed to trace rays across a 2D geometry. Each surface
 *          type to be used in a simulation is a subclass of the Surface
 *          class
 */
class Surface {
protected:
    /** The surface's name */
    char* _surface_name;
    /** A static class variable to generate a UID for each new surface */
    static int _n;
    /** The surface's unique identifier */
    int _uid;
    /** The type of surface */
    surfaceType _surface_type;
    /** The boundary condition for this Surface */
    boundaryType _boundary_type;
public:
    Surface(const char* surface_name=(char*)"");
    virtual ~Surface();

    char* getSurfaceName();
    int getUid() const;    
    boundaryType getBoundaryType() const;
    surfaceType getSurfaceType() const;

    void setBoundaryType(boundaryType type);

    /**
     * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y,z) \f$ 
     *        with respect to a quadratic surface \f$ f(x,y,z) \f$.
     * @param neutron the neutron of interest
     *
     */
    virtual float evaluate(neutron* neutron) =0;

    /**
     * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y,z) \f$ 
     *        with respect to a quadratic surface \f$ f(x,y,z) \f$.
     * @param x the x-coordinate of interest
     * @param y the y-coordinate of interest
     * @param z the z-coordinate of interest
     *
     */
    virtual float evaluate(float x, float y, float z) =0;

    /**
     * @brief Computes the nearest distance between a neutron at some
     *        location and this surface.
     * @details This virtual class method must be implemented for each 
     *          surface type to be used in a PINSPEC simulation.
     * @param neutron the neutron of interest
     */ 
    virtual float computeNearestDistance(neutron* neutron) =0;

    /**
     * @brief Determines whether or not a neutron at some location is
     *        on the surface.
     * @details This virtual class method must be implemented for each
     *          surface type to be used in a PINSPEC simulation.
     * @param neutron the neutron of interest
     */
    virtual bool onSurface(neutron* neutron) =0;

    /**
     * @brief Perfectly reflects a neutron at the surface for REFLECTIVE
     *        boundary conditions.
     * @details This virtual class method must be implemented for each
     *          surface type to be used in a PINSPEC simulation.
     * @param neutron the neutron of interest
     */
    virtual void reflectNeutron(neutron* neutron) =0;
};


/**
 * @class XPlane Surface.h "pinspec/src/Surface.h"
 * @brief The XPlane is a plane perpendicular to the y-axis.
 * @details The XPlane represents planes perpendicular to the x-axis using the 
 *          quadratic surface formulation. The XPlane is a Surface subclass.
 */
class XPlane: public Surface {
private:
    /** The location of the plane's intersection with the x-axis */
    float _x;
public:
    XPlane(const char* surface_name=(char*)"");
    virtual ~XPlane();
    float getX();
    void setX(float x);
    float evaluate(neutron* neutron);
    float evaluate(float x, float y, float z);
    float computeNearestDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    void reflectNeutron(neutron* neutron);
};


/**
 * @class YPlane Surface.h "pinspec/src/Surface.h"
 * @brief The YPlane is a a plane perpendicular to the y-axis.
 * @details The YPlane represents planes perpendicular to the y-axis using the
 *          quadratic surface formulation. The YPlane is a Surface subclass.
 */ 
class YPlane: public Surface {
private:
    /** The location of the plane's intersection with the y-axis */
    float _y;
public:
    YPlane(const char* surface_name=(char*)"");
    virtual ~YPlane();
    float getY();
    void setY(float y);
    float evaluate(neutron* neutron);
    float evaluate(float x, float y, float z);
    float computeNearestDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    void reflectNeutron(neutron* neutron);
};


/**
 * @class ZCylinder Surface.h "pinspec/src/Surface.h"
 * @brief The ZCylinder is the locus of a point equidistant from a fixed point.
 * @details The ZCylinder represents a set of points within a given distance
 *          from a point in the xy-plane called the cylinder's center.
 */
class ZCylinder: public Surface {
protected:
    /** The cylinder's radius */
    float _r;
    /** The square of the cylinder's radius */
    float  _r_squared;
    /** The x-coordinate of the cylinder's center */
    float _x0;
    /** The y-coordinate of the cylinder's center */
    float _y0;
public:
    ZCylinder(const char* surface_name=(char*)"");
    virtual ~ZCylinder();
    float getX0();
    float getY0();
    float getRadius();
    void setX0(float x0);
    void setY0(float y0);
    void setRadius(float r);
    float evaluate(neutron* neutron);
    float evaluate(float x, float y, float z);
    float computeNearestDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    void reflectNeutron(neutron* neutron);
};

#endif /* SURFACE_H_ */
