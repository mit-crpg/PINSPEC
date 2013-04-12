/**
 * @file RegionFactory.h
 * @brief The RegionFactory class.
 * @author William Boyd (wboyd@mit.edu)
 * @date April 12, 2013
 */

#ifndef REGIONFACTORY_H_
#define REGIONFACTORY_H_

#include "Region.h"


/**
 * @class RegionFactory RegionFactory.h "pinspec/src/RegionFactory.h"
 * @brief A utility class for creating instances of regions.
 */
class RegionFactory {

private:
    /**
     * @brief RegionFactory constructor.
     */
    RegionFactory() { }

    /**
     * @brief Assignment operator for static referencing of the RegionFactory.
     * @param & the RegionFactory static class object
     * @return a pointer to the RegionFactory static class object
     */
    RegionFactory &operator=(const RegionFactory &) { return *this; }
 
public:
    ~RegionFactory() { }

    /**
     * @brief Returns a static instance of the RegionFactory class.
     * @return a pointer to the static RegionFactory class
     */
    static RegionFactory *Get() {
        static RegionFactory instance;
	return &instance;
    }

    Region* createRegion(regionType region_type, 
		       const char* region_name=(char*)"");
};

#endif /* REGIONFACTORY_H_ */
