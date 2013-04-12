#include "RegionFactory.h"

/**
 * @brief Method to create a region for a PINSPEC monte carlo simulation.
 * @param region_type the type of region (ie, INFINITE_HOMOGENEOUS, 
 *        HOMOGENEOUS_EQUIVALENCE, or HETEROGENEOUS)
 * @param region_name a character array for the name of the region
 * @return a pointer to the newly created region
 */
Region* RegionFactory::createRegion(regionType region_type, 
				    const char* region_name) {
    if (region_type == INFINITE_MEDIUM)
        return new InfiniteMediumRegion(region_name);
    else if (region_type == EQUIVALENT_FUEL)
        return new EquivalenceRegion(region_type, region_name);
    else if (region_type == EQUIVALENT_MODERATOR)
        return new EquivalenceRegion(region_type, region_name);
    else if (region_type == BOUNDED_FUEL)
        return new BoundedFuelRegion(region_name);
    else if (region_type == BOUNDED_MODERATOR)
        return new BoundedModeratorRegion(region_name);
    else
        return new BoundedRegion(region_name);
}
