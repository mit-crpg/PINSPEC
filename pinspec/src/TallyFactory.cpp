/*
 * TallyFactory.cpp
 *
 *  Created on: Mar 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "TallyFactory.h"


Tally* TallyFactory::createTally(char* tally_name, Isotope* isotope, 
														tallyType tally_type) {
	if (tally_type == FLUX)
		log_printf(ERROR, "Unable to create a FLUX type Tally for an "
				"Isotope. FLUX Tallies are only supported for Regions.");
	if (tally_type == LEAKAGE_RATE)
		log_printf(ERROR, "Unable to create a LEAKAGE_RATE type tally for an"
			" Isotope. LEAKAGE_RATE Tallies are only supported for Regions.");

	if (tally_type == COLLISION_RATE)
		return new IsotopeCollisionRateTally(tally_name, isotope);
	else if (tally_type == ELASTIC_RATE)
		return new IsotopeElasticRateTally(tally_name, isotope);
	else if (tally_type == ABSORPTION_RATE)
		return new IsotopeAbsorptionRateTally(tally_name, isotope);
	else if (tally_type == CAPTURE_RATE)
		return new IsotopeCaptureRateTally(tally_name, isotope);
	else if (tally_type == FISSION_RATE)
		return new IsotopeFissionRateTally(tally_name, isotope);
	else if (tally_type == TRANSPORT_RATE)
		return new IsotopeTransportRateTally(tally_name, isotope);
	else
		return new IsotopeDiffusionRateTally(tally_name, isotope);
}


Tally* TallyFactory::createTally(char* tally_name, Material* material, 
														tallyType tally_type) {


	if (tally_type == FLUX)
		return new MaterialFluxTally(tally_name, material);
	else if (tally_type == LEAKAGE_RATE)
		return new MaterialLeakageRateTally(tally_name, material);
	else if (tally_type == COLLISION_RATE)
		return new MaterialCollisionRateTally(tally_name, material);
	else if (tally_type == ELASTIC_RATE)
		return new MaterialElasticRateTally(tally_name, material);
	else if (tally_type == ABSORPTION_RATE)
		return new MaterialAbsorptionRateTally(tally_name, material);
	else if (tally_type == CAPTURE_RATE)
		return new MaterialCaptureRateTally(tally_name, material);
	else if (tally_type == FISSION_RATE)
		return new MaterialFissionRateTally(tally_name, material);
	else if (tally_type == TRANSPORT_RATE)
		return new MaterialTransportRateTally(tally_name, material);
	else
		return new MaterialDiffusionRateTally(tally_name, material);
}


Tally* TallyFactory::createTally(char* tally_name, Region* region, 
														tallyType tally_type) {

	if (tally_type == FLUX)
		return new RegionFluxTally(tally_name, region);
	else if (tally_type == LEAKAGE_RATE)
		return new RegionLeakageRateTally(tally_name, region);
	else if (tally_type == COLLISION_RATE)
		return new RegionCollisionRateTally(tally_name, region);
	else if (tally_type == ELASTIC_RATE)
		return new RegionElasticRateTally(tally_name, region);
	else if (tally_type == ABSORPTION_RATE)
		return new RegionAbsorptionRateTally(tally_name, region);
	else if (tally_type == CAPTURE_RATE)
		return new RegionCaptureRateTally(tally_name, region);
	else if (tally_type == FISSION_RATE)
		return new RegionFissionRateTally(tally_name, region);
	else if (tally_type == TRANSPORT_RATE)
		return new RegionTransportRateTally(tally_name, region);
	else
		return new RegionDiffusionRateTally(tally_name, region);
}


Tally* TallyFactory::createTally(char* tally_name, Geometry* geometry, 
														tallyType tally_type) {

	if (tally_type == FLUX)
		return new GeometryFluxTally(tally_name, geometry);
	else if (tally_type == LEAKAGE_RATE)
		return new GeometryLeakageRateTally(tally_name, geometry);
	else if (tally_type == COLLISION_RATE)
		return new GeometryCollisionRateTally(tally_name, geometry);
	else if (tally_type == ELASTIC_RATE)
		return new GeometryElasticRateTally(tally_name, geometry);
	else if (tally_type == ABSORPTION_RATE)
		return new GeometryAbsorptionRateTally(tally_name, geometry);
	else if (tally_type == CAPTURE_RATE)
		return new GeometryCaptureRateTally(tally_name, geometry);
	else if (tally_type == FISSION_RATE)
		return new GeometryFissionRateTally(tally_name, geometry);
	else if (tally_type == TRANSPORT_RATE)
		return new GeometryTransportRateTally(tally_name, geometry);
	else
		return new GeometryDiffusionRateTally(tally_name, geometry);

}
