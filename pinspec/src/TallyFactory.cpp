#include "TallyFactory.h"

/**
 * @brief Method to create a tally for some tally type within an 
 *        isotope.
 * @param isotope a pointer to the isotope to tally
 * @param tally_type the type of tally (ie, FLUX, CAPTURE_RATE, etc.)
 * @param tally_name a character array for the name of the tally
 * @return a pointer to the newly created tally
 */
Tally* TallyFactory::createTally(Isotope* isotope, tallyType tally_type, 
                                                  char* tally_name) {
    if (tally_type == FLUX)
        log_printf(ERROR, "Unable to create a FLUX type Tally for an "
		   "Isotope. FLUX tallies are only supported for "
		   "materials, regions and the geometry.");
    if (tally_type == LEAKAGE_RATE)
        log_printf(ERROR, "Unable to create a LEAKAGE_RATE type tally for an"
		   " Isotope. LEAKAGE_RATE tallies are only supported for "
		   "materials, regions and the geometry.");
    if (tally_type == INTERCOLLISION_TIME)
        log_printf(ERROR, "Unable to create an INTERCOLLISION_TIME type "
		   "tally for an Isotope. INTERCOLLISION_TIME tallies are only "
		   "supported for materials, regions and the geometry.");
    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
                                                        " TallyFactory.");

    if (tally_type == COLLISION_RATE)
        return new IsotopeCollisionRateTally(isotope, tally_name);
    else if (tally_type == ELASTIC_RATE)
        return new IsotopeElasticRateTally(isotope, tally_name);
    else if (tally_type == GROUP_TO_GROUP_RATE)
        return new IsotopeGroupRateTally(isotope, tally_name);
    else if (tally_type == OUTSCATTER_RATE)
        return new IsotopeOutScatterRateTally(isotope, tally_name);
    else if (tally_type == ABSORPTION_RATE)
        return new IsotopeAbsorptionRateTally(isotope, tally_name);
    else if (tally_type == CAPTURE_RATE)
        return new IsotopeCaptureRateTally(isotope, tally_name);
    else if (tally_type == FISSION_RATE)
        return new IsotopeFissionRateTally(isotope, tally_name);
    else if (tally_type == TRANSPORT_RATE)
        return new IsotopeTransportRateTally(isotope, tally_name);
    else
        return new IsotopeDiffusionRateTally(isotope, tally_name);
}


/**
 * @brief Method to create a tally for some tally type within a material.
 * @param material a pointer to the material within which to tally
 * @param tally_type the type of tally (ie, FLUX, CAPTURE_RATE, etc.)
 * @param tally_name a character array for the name of the tally
 * @return a pointer to the newly created tally
 */
Tally* TallyFactory::createTally(Material* material, tallyType tally_type,
                                                  char* tally_name) {

    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
                                                        " TallyFactory.");

    if (tally_type == FLUX)
        return new MaterialFluxTally(material, tally_name);
    else if (tally_type == LEAKAGE_RATE)
        return new MaterialLeakageRateTally(material, tally_name);
    else if (tally_type == INTERCOLLISION_TIME)
        return new MaterialInterCollisionTimeTally(material, tally_name);    
    else if (tally_type == COLLISION_RATE)
        return new MaterialCollisionRateTally(material, tally_name);
    else if (tally_type == ELASTIC_RATE)
        return new MaterialElasticRateTally(material, tally_name);
    else if (tally_type == GROUP_TO_GROUP_RATE)
        return new MaterialGroupRateTally(material, tally_name);
    else if (tally_type == OUTSCATTER_RATE)
        return new MaterialOutScatterRateTally(material, tally_name);
    else if (tally_type == ABSORPTION_RATE)
        return new MaterialAbsorptionRateTally(material, tally_name);
    else if (tally_type == CAPTURE_RATE)
        return new MaterialCaptureRateTally(material, tally_name);
    else if (tally_type == FISSION_RATE)
        return new MaterialFissionRateTally(material, tally_name);
    else if (tally_type == TRANSPORT_RATE)
        return new MaterialTransportRateTally(material, tally_name);
    else
        return new MaterialDiffusionRateTally(material, tally_name);
}


/**
 * @brief Method to create a tally for some tally type within a region.
 * @param region a pointer to the region within which to tally
 * @param tally_type the type of tally (ie, FLUX, CAPTURE_RATE, etc.)
 * @param tally_name a character array for the name of the tally
 * @return a pointer to the newly created tally
 */
Tally* TallyFactory::createTally(Region* region, tallyType tally_type,
                                                  char* tally_name) {

    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
		   " TallyFactory.");

    if (tally_type == FLUX)
        return new RegionFluxTally(region, tally_name);
    else if (tally_type == LEAKAGE_RATE)
        return new RegionLeakageRateTally(region, tally_name);
    else if (tally_type == INTERCOLLISION_TIME)
        return new RegionInterCollisionTimeTally(region, tally_name);    
    else if (tally_type == COLLISION_RATE)
        return new RegionCollisionRateTally(region, tally_name);
    else if (tally_type == ELASTIC_RATE)
        return new RegionElasticRateTally(region, tally_name);
    else if (tally_type == GROUP_TO_GROUP_RATE)
	return new RegionGroupRateTally(region, tally_name);
    else if (tally_type == OUTSCATTER_RATE)
	return new RegionOutScatterRateTally(region, tally_name);
    else if (tally_type == ABSORPTION_RATE)
        return new RegionAbsorptionRateTally(region, tally_name);
    else if (tally_type == CAPTURE_RATE)
        return new RegionCaptureRateTally(region, tally_name);
    else if (tally_type == FISSION_RATE)
        return new RegionFissionRateTally(region, tally_name);
    else if (tally_type == TRANSPORT_RATE)
        return new RegionTransportRateTally(region, tally_name);
    else
        return new RegionDiffusionRateTally(region, tally_name);
}


/**
 * @brief Method to create a tally for some tally type within the geometry.
 * @param geometry a pointer to the geometry within which to tally
 * @param tally_type the type of tally (ie, FLUX, CAPTURE_RATE, etc.)
 * @param tally_name a character array for the name of the tally
 * @return a pointer to the newly created tally
 */
Tally* TallyFactory::createTally(Geometry* geometry, tallyType tally_type,
                                                  char* tally_name) {

    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
                                                        " TallyFactory.");

    if (tally_type == FLUX)
        return new GeometryFluxTally(geometry, tally_name);
    else if (tally_type == LEAKAGE_RATE)
        return new GeometryLeakageRateTally(geometry, tally_name);
    else if (tally_type == INTERCOLLISION_TIME)
        return new GeometryInterCollisionTimeTally(geometry, tally_name);    
    else if (tally_type == COLLISION_RATE)
        return new GeometryCollisionRateTally(geometry, tally_name);
    else if (tally_type == ELASTIC_RATE)
        return new GeometryElasticRateTally(geometry, tally_name);
    else if (tally_type == GROUP_TO_GROUP_RATE)
	return new GeometryGroupRateTally(geometry, tally_name);
    else if (tally_type == OUTSCATTER_RATE)
	return new GeometryOutScatterRateTally(geometry, tally_name);
    else if (tally_type == ABSORPTION_RATE)
        return new GeometryAbsorptionRateTally(geometry, tally_name);
    else if (tally_type == CAPTURE_RATE)
        return new GeometryCaptureRateTally(geometry, tally_name);
    else if (tally_type == FISSION_RATE)
        return new GeometryFissionRateTally(geometry, tally_name);
    else if (tally_type == TRANSPORT_RATE)
        return new GeometryTransportRateTally(geometry, tally_name);
    else
        return new GeometryDiffusionRateTally(geometry, tally_name);
}
