/*
 * TallyFactory.cpp
 *
 *  Created on: Mar 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "TallyFactory.h"


Tally* TallyFactory::createTally(Isotope* isotope, tallyType tally_type, 
                                                  const char* tally_name) {
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
		return new IsotopeElasticRateTally(tally_name, isotope);
	else if (tally_type == OUTSCATTER_RATE)
	    return new IsotopeOutscatterRateTally(tally_name, isotope);
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


Tally* TallyFactory::createTally(Material* material, tallyType tally_type,
                                                  const char* tally_name) {

    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
                                                        " TallyFactory.");

	if (tally_type == FLUX)
		return new MaterialFluxTally(material, tally_name);
	else if (tally_type == LEAKAGE_RATE)
		return new MaterialLeakageRateTally(tally_name, material);
	else if (tally_type == INTERCOLLISION_TIME)
		return new MaterialInterCollisionTimeTally(tally_name, material);    
	else if (tally_type == COLLISION_RATE)
		return new MaterialCollisionRateTally(material, tally_name);
	else if (tally_type == ELASTIC_RATE)
		return new MaterialElasticRateTally(tally_name, material);
	else if (tally_type == OUTSCATTER_RATE)
	    return new MaterialOutscatterRateTally(tally_name, material);
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


Tally* TallyFactory::createTally(Region* region, tallyType tally_type,
                                                  const char* tally_name) {

    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
                                                        " TallyFactory.");

	if (tally_type == FLUX)
		return new RegionFluxTally(region, tally_name);
	else if (tally_type == LEAKAGE_RATE)
		return new RegionLeakageRateTally(region, tally_name);
    else if (tally_type == INTERCOLLISION_TIME)
		return new RegionInterCollisionTimeTally(tally_name, region);    	else if (tally_type == OUTSCATTER_RATE)
	    return new RegionOutscatterRateTally(tally_name, region); 
	else if (tally_type == COLLISION_RATE)
		return new RegionCollisionRateTally(region, tally_name);
	else if (tally_type == ELASTIC_RATE)
		return new RegionElasticRateTally(region, tally_name);
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


Tally* TallyFactory::createTally(Geometry* geometry, tallyType tally_type,
                                                  const char* tally_name) {

    if (tally_type == DERIVED)
        log_printf(ERROR, "DERIVED type tallies cannot be created by the "
                                                        " TallyFactory.");

	if (tally_type == FLUX)
		return new GeometryFluxTally(geometry, tally_name);
	else if (tally_type == LEAKAGE_RATE)
		return new GeometryLeakageRateTally(geometry, tally_name);
    else if (tally_type == INTERCOLLISION_TIME)
		return new GeometryInterCollisionTimeTally(tally_name, geometry);    	
    else if (tally_type == OUTSCATTER_RATE)
	    return new GeometryOutscatterRateTally(tally_name, geometry);
	else if (tally_type == COLLISION_RATE)
		return new GeometryCollisionRateTally(geometry, tally_name);
	else if (tally_type == ELASTIC_RATE)
		return new GeometryElasticRateTally(geometry, tally_name);
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
