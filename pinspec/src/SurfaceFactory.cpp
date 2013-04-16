#include "SurfaceFactory.h"

/**
 * @brief Method to create a bounding surface for a region.
 * @param surface_type the type of surface (ie, XPLANE, YPLANE, ZCYLINDER, etc.)
 * @param surface_name a character array for the name of the surface
 * @return a pointer to the newly created surface
 */
Surface* SurfaceFactory::createSurface(surfaceType surface_type,                                                   const char* surface_name) {
    if (surface_type == XPLANE)
        return new XPlane(surface_name);
    else if (surface_type == YPLANE)
        return new YPlane(surface_name);
    else
        return new ZCylinder(surface_name);
}
