#include "Surface.h"

int Surface::_n = 1;

/**
 * @brief Surface class constructor.
 * @param surface_name the (optional) name of the surface
 * @details By default, the constructor sets the surface boundary conditions
 *          to vacuum.
 */
Surface::Surface(const char* surface_name) {
    _surface_name = (char*)surface_name;
    _uid = _n;
    _n++;
    _boundary_type = VACUUM;
}


/**
 * @brief Surface destructor.
 */
Surface::~Surface() { }


/**
 * @brief Return the name of the surface.
 * @return a character array representing this surface's name
 */
char* Surface::getSurfaceName() {
    return _surface_name;
}


/**
 * @brief Returns the unique ID auto-generated for the surface.
 * @return a unique ID for the surface
 */
int Surface::getUid() const {
    return _uid;
}


/**
 * @brief Returns this surface's type (XPLANE, YPLANE, CIRCLE, etc.)
 * @returns Surface type
 */
surfaceType Surface::getSurfaceType() const {
    return _surface_type;
}


/**
 * @brief Returns this surface's boundary type.
 * @details Returns the surface's boundary type which can be REFLECTIVE, VACUUM, 
 *          or INTERFACE.
 * @return Surface boundary type
 */
boundaryType Surface::getBoundaryType() const {
    return _boundary_type;
}


/**
 * @brief Sets the boundary type for this Surface. 
 * @details Sets the Surface's boundary type which can be REFLECTIVE, VACUUM, 
 *          or INTERFACE).
 * @param type the boundary type
 */
void Surface::setBoundaryType(boundaryType type) {
    _boundary_type = type;
}


/******************************************************************************
 *****************************   XPlane   *************************************
 *****************************************************************************/

/**
 * @brief The XPlane constructor.
 * @param surface_name the (optional) name of the surface
 * @details Assigns a default value of x=0 for the x-axis intersection point.
 */
XPlane::XPlane(const char* surface_name) : Surface(surface_name) {
    _x = 0.0;
};


/**
 * @brief XPlane destructor.
 */
XPlane::~XPlane() { };


/**
 * @brief Returns the x-coordinate of the plane's intersection point with the 
 *        x-axis.
 * @return the x-coordinate of the intersection point
 */
float XPlane::getX() {
    return _x;
}


/**
 * @brief Sets the x-coordinate of the plane's intersection point with the 
 *        x-axis.
 * @param x the x-coordinate of the intersectiont point
 */
void XPlane::setX(float x) {
    _x = x;
}


/**
 * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y,z) \f$ 
 *        with respect to the quadratic surface representing this x-plane:
 *        \f$ f(x,y,z) = qx + d = x - x_0 \f$
 * @param neutron the neutron of interest
 */
float XPlane::evaluate(neutron* neutron) {
    float x = neutron->_x;
    return (x - _x);
}


/**
 * @brief Returns the evaluation of a 3D coordinate \f$ (x,y,z) \f$ 
 *        with respect to the quadratic surface representing this x-plane:
 *        \f$ f(x,y,z) = qx + d = x - x_0 \f$
 * @param x the x-coordinate of interest
 * @param y the y-coordinate of interest
 * @param z the z-coordinate of interest
 */
float XPlane::evaluate(float x, float y, float z) {
    return (x - _x);
}



/**
 * @brief Computes the parametrized distance to the XPlane along a neutron's 
 *        trajectory.
 * @details The distance returned is a parametrized distance along the unit
 *          vectory defining the neutron's trajectory, not a Cartesian
 *          distance. If the neutron is traveling away from the surface then
 *          the distance returned will be infinity. The equation used to 
 *          compute the parametrized distance \f$ d \f$ from a neutron at 
 *          \f$ x \f$ traveling with component an x-component of its unit 
 *          trajectory vector \f$ u \f$ to this xplane at \f$ x_0 \f$ is:
 *          \f$ d = \frac{x_0 - x}{u} \f$
 *          
 * @param neutron a pointer to a neutron struct
 * @return the parametrized distance to the surface along the neutron's trajectory
 */
float XPlane::computeParametrizedDistance(neutron* neutron) {

    /* Set dist to infinity to begin with */
    float dist = std::numeric_limits<int>::max();
    float x = neutron->_x;
    float u = neutron->_u;

    /* If the neutron is traveling towards the surface, return the distance */
    if ((x - _x < -1E-6 && u > 0.0) || (x - _x > 1E-6 && u < 0.0))
        dist = (_x - x) / u;

    /* If the neutron is not traveling towards the surface, return infinity */
    return dist;
}


/**
 * @brief Checks whether a neutron is on the XPlane.
 * @details The threshold used to compute whether or not a neutron is on the
 *          on the XPlane is 1E-6 for the difference between the x-coordinate 
 *          of the neutron's position and the location of the XPlane.
 * @param neutron the neutron of interest
 * @return true if on the XPlane, otherwise false
 */
bool XPlane::onSurface(neutron* neutron) {
    if (fabs(_x - neutron->_x) < 1E-6)
        return true;

    return false;
}


/**
 * @brief Perfectly reflects a neutron at a xplane.
 * @param neutron the neutron of interest
 */
void XPlane::reflectNeutron(neutron* neutron) {
    /* Reverse the x-component of the velocity vector */
    neutron->_u *= -1.0;
}



/******************************************************************************
 *****************************   YPlane   *************************************
 *****************************************************************************/

/**
 * @brief The YPlane constructor.
 * @param surface_name the (optional) name of the surface
 * @details Assigns a default value of y=0 for the y-axis intersection point.
 */
YPlane::YPlane(const char* surface_name) : Surface(surface_name) {
    _y = 0.0;
};


/**
 * @brief YPlane destructor.
 */
YPlane::~YPlane() { };



/**
 * @brief Returns the y-coordinate of the plane's intersection point with the 
 *        y-axis.
 * @return the y-coordinate of the intersection point
 */
float YPlane::getY() {
    return _y;
}


/**
 * @brief Sets the x-coordinate of the plane's intersection point with the 
 *        y-axis.
 * @param y the y-coordinate of the intersectiont point
 */
void YPlane::setY(float y) {
    _y = y;
}


/**
 * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y) \f$ 
 *        with respect to the quadratic surface representing this y-plane:
 *        \f$ f(x,y) = qy + d = y - y_0 \f$
 * @param neutron the neutron of interest
 */
float YPlane::evaluate(neutron* neutron) {
    float y = neutron->_y;
    return (y - _y);
}


/**
 * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y,z) \f$ 
 *        with respect to the quadratic surface representing this y-plane:
 *        \f$ f(x,y,z) = qy + d = y - y_0 \f$
 * @param x the x-coordinate of interest
 * @param y the y-coordinate of interest
 * @param z the z-coordinate of interest 
 */
float YPlane::evaluate(float x, float y, float z) {
    return (y - _y);
}


/**
 * @brief Computes the parametrized distance to the YPlane along a neutron's 
 *        trajectory.
 * @details The distance returned is a parametrized distance along the unit
 *          vectory defining the neutron's trajectory, not a Cartesian
 *          distance. If the neutron is traveling away from the surface then
 *          the distance returned will be infinity. The equation used to 
 *          compute the parametrized distance \f$ d \f$ from a neutron at 
 *          \f$ y \f$ traveling with component an y-component of its unit 
 *          trajectory vector \f$ u \f$ to this xplane at \f$ y_0 \f$ is:
 *          \f$ d = \frac{y_0 - y}{u} \f$
 *          
 * @param neutron a pointer to a neutron struct
 * @return the parametrized distance to the surface along the neutron's trajectory
 */
float YPlane::computeParametrizedDistance(neutron* neutron) {

    /* Set dist to infinity to begin with */
    float dist = std::numeric_limits<int>::max();
    float y = neutron->_y;
    float v = neutron->_v;

    /* If the neutron is traveling towards the surface, return the distance */
    if ((y - _y < -1E-6 && v > 0.0) || (y - _y > 1E-6 && v < 0.0))
        dist = (_y - y) / v;

    /* If the neutron is not traveling towards the surface, return infinity */
    return dist;
}


/**
 * @brief Checks whether a neutron is on the YPlane.
 * @details The threshold used to compute whether or not a neutron is on the
 *          on the YPlane is 1E-6 for the difference between the y-coordinate 
 *          of the neutron's position and the location of the YPlane
 * @param neutron the neutron of interest to check
 * @return true if on the YPlane, otherwise false
 */
bool YPlane::onSurface(neutron* neutron) {
    if (fabs(_y - neutron->_y) < 1E-6)
        return true;

    return false;
}



/**
 * @brief Perfectly reflects a neutron at a yplane.
 * @param neutron the neutron of interest
 */
void YPlane::reflectNeutron(neutron* neutron) {
    neutron->_v *= -1.0;
}


/******************************************************************************
 ****************************   ZCylinder   ***********************************
 *****************************************************************************/

/**
 * @brief The ZCylinder constructor.
 * @param surface_name the (optional) name of the surface
 * @details Assigns default values for the center of the cylinder (x=0, y=0)
 *          and a radius of 0.
 */
ZCylinder::ZCylinder(const char* surface_name) : Surface(surface_name) {
    _x0 = 0;
    _y0 = 0;
    _r = 0;
    _r_squared = 0;
};


/**
 * @brief ZCylinder destructor.
 */
ZCylinder::~ZCylinder() { };


/**
 * @brief Returns the x-coordinate of the cylinder's center.
 * @return the x-coordinate of the cylinder center
 */
float ZCylinder::getX0() {
    return _x0;
}


/**
 * @brief Returns the y-coordinate of the cylinder's center.
 * @return the y-coordinate of the cylinder center
 */
float ZCylinder::getY0() {
    return _y0;
}


/**
 * @brief Returns the radius of the cylinder.
 * @return the cylinder radius
 */
float ZCylinder::getRadius() {
    return _r;
}


/**
 * @brief Sets the x-coordinate of the cylinder's center.
 * @param x0 the x-coordinate of the cylinder center
 */
void ZCylinder::setX0(float x0) {
     _x0 = x0;
}


/**
 * @brief Sets the y-coordinate of the cylinder's center.
 * @param y0 the y-coordinate of the cylinder center
 */
void ZCylinder::setY0(float y0) {
    _y0 = y0;
}


/**
 * @brief Sets the radius of the cylinder's center.
 * @param r the cylinder's radius
 */
void ZCylinder::setRadius(float r) {
     _r = r;
     _r_squared = r*r;
}


/**
 * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y) \f$ 
 *        with respect to the quadratic surface representing this cylinder:
 *        \f$ f(x,y) = (x - x_0)^2 + (y - y_0)^2 - R^2 \f$
 * @param neutron the neutron of interest
 */
float ZCylinder::evaluate(neutron* neutron) {
    float x = neutron->_x;
    float y = neutron->_y;
    return (x - _x0) * (x - _x0) + (y - _y0) * (y - _y0) - _r_squared;
}


/**
 * @brief Returns the evaluation of a neutron's coordinates \f$ (x,y,z) \f$ 
 *        with respect to the quadratic surface representing this cylinder:
 *        \f$ f(x,y,z) = (x - x_0)^2 + (y - y_0)^2 - R^2 \f$
 * @param x the x-coordinate of interest
 * @param y the y-coordinate of interest
 * @param z the z-coordinate of interest
 */
float ZCylinder::evaluate(float x, float y, float z) {
    return (x - _x0) * (x - _x0) + (y - _y0) * (y - _y0) - _r_squared;
}


/**
 * @brief Computes the parametrized distance to the ZCylinder along a neutron's 
 *        trajectory.
 * @details The distance returned is a parametrized distance along the unit
 *          vectory defining the neutron's trajectory, not a Cartesian
 *          distance. If the neutron is traveling away from the surface then
 *          the distance returned will be infinity. This method solves a
 *          quadratic equation for 0, 1, or 2 surface intersection points
 *          and selects the one nearest to the neutron.
 *          
 * @param neutron a pointer to a neutron struct
 * @return the parametrized distance to the surface along the neutron's trajectory
 */
float ZCylinder::computeParametrizedDistance(neutron* neutron) {

    float x = neutron->_x;
    float y = neutron->_y;
    float u = neutron->_u;
    float v = neutron->_v;

    float r_squared = x*x + y*y;

    /* Compute temporary variables for each term in the quadratic equation */
    float a = u*u + v*v;
    float b = 2.0*x*u - 2.0*_x0*u + 2.0*y*v - 2.0*_y0*v;
    float c = r_squared + _x0*_x0 - 2.0*_x0*x + _y0*_y0 - 2.0*_y0*y 
                  - _r_squared;

    /* Compute the discriminant in the quadratic formula */
    float discr = b*b - 4.0*a*c;

    /* There is not an intersection point */
    if (discr < 0.0)
        return std::numeric_limits<int>::max();

    /* There is one intersection point */
    else if (discr == 0.0) {

        float dist = -b / (2.0*a);

        if (dist > 0.0)
            return dist;
        else
            return std::numeric_limits<int>::max();
    }

    /* There are two intersection points */
    else {
        discr = sqrt(discr);
        float dist1 = (-b + discr) / (2.0*a);
        float dist2 = (-b - discr) / (2.0*a);

        if (dist1 < 0.0 && dist2 < 0.0)
            return std::numeric_limits<int>::max();

        /* Determine which distance is shorter */
        if (dist1 > 0.0 && dist2 < 0.0)
	    return dist1;
	else if (dist2 > 0.0 && dist1 < 0.0)
	    return dist2;
        else if (dist1 < dist2 && dist1 > 0.0)
	    return dist1;
	else
	    return dist2;
    }
}


/**
 * @brief Checks whether a neutron is on the ZCylinder.
 * @details The threshold used to compute whether or not a neutron is on the
 *          on the neutron is 1E-6 for the difference between the distance
 *          between the neutron and the cylinder center and the radius of the
 *          cylinder.
 * @param neutron the neutron of interest
 * @return true if on the XPlane, otherwise false
 */
bool ZCylinder::onSurface(neutron* neutron) {

    float r_squared = (neutron->_y - _y0) * (neutron->_y - _y0) + 
                       (neutron->_x - _x0) * (neutron->_x - _x0);

    if (fabs(_r_squared - r_squared) < 1E-6)
        return true;
    else
        return false;
}


/**
 * @brief Perfectly reflects a neutron at a cylinder.
 * @param neutron the neutron of interest
 */
void ZCylinder::reflectNeutron(neutron* neutron) {

    /* Compute the vector normal to the cylinder at the intersection point of the
     * neutron's trajectory and the cylinder */
    float x1 = neutron->_x - _x0;
    float y1 = neutron->_y - _y0;
   
    /* Compute the unit vector of the neutron's trajectory */
    float u = neutron->_u;
    float v = neutron->_v;

    /* Compute the angle between the two vectors at the reflection point on 
     *the cylinder's surface */
    float theta = acos(dotProduct2D(x1, y1, u, v) 
			/ (norm2D(x1, y1) * norm2D(u, v)));

    /* Reflect the particle around the vector normal to the cylinder surface */
    float rotation_angle = M_PI - 2.0 * theta;

    /* Update the neutron's trajectory vector using the 2D rotation matrix */
    neutron->_u = neutron->_u * cos(rotation_angle) -
                  neutron->_v * sin(rotation_angle);
    neutron->_v = neutron->_u * sin(rotation_angle) +
                  neutron->_v * cos(rotation_angle);

    return;
}
