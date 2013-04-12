/**
 * @file SurfaceFactory.h
 * @brief The SurfaceFactory class.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 19, 2013
 */

#ifndef SURFACEFACTORY_H_
#define SURFACEFACTORY_H_

#include "Surface.h"


/**
 * @class SurfaceFactory SurfaceFactory.h "pinspec/src/SurfaceFactory.h"
 * @brief A utility class for creating instances of surface.
 */
class SurfaceFactory {

private:
    /**
     * @brief SurfaceFactory constructor.
     */
    SurfaceFactory() { }

    /**
     * @brief Assignment operator for static referencing of the SurfaceFactory.
     * @param & the SurfaceFactory static class object
     * @return a pointer to the SurfaceFactory static class object
     */
    SurfaceFactory &operator=(const SurfaceFactory &) { return *this; }
 
public:
    ~SurfaceFactory() { }

    /**
     * @brief Returns a static instance of the SurfaceFactory class.
     * @return a pointer to the static SurfaceFactory class
     */
    static SurfaceFactory *Get() {
        static SurfaceFactory instance;
	return &instance;
    }

    Surface* createSurface(surfaceType surface_type, 
			   const char* surface_name=(char*)""); 
};

#endif /* SURFACEFACTORY_H_ */
