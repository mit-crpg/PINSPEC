/*
 * Tally.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TALLY_H_
#define TALLY_H_

#ifdef __cplusplus
#include <limits>
#include <math.h>
#include "log.h"
#include "arraycreator.h"
#include "Neutron.h"
#include "Material.h"
#include "Isotope.h"
#include "Region.h"

class Geometry;


#define NEUTRON_MASS 939565378           /* Mass of neutron in eV / c^2 */
#define LIGHT_SPEED 299792458               /* Speed of light in m / s */


/* The domain in which this tally resides */
typedef enum tallyDomainTypes {
	MATERIAL,
	ISOTOPE,
	REGION,
	GEOMETRY
} tallyDomainType;


/* The type of precision trigger for a Tally object - represents
 * what kind of precision we should track to determine when to end
 * the simulation */
typedef enum triggerTypes {
    VARIANCE,
    STANDARD_DEVIATION,
    RELATIVE_ERROR,
    NONE
} triggerType;


/* Type of tallies */
typedef enum tallyTypes {
	FLUX,
	LEAKAGE_RATE,
	COLLISION_RATE,
	INTERCOLLISION_TIME,
	ELASTIC_RATE,
	ABSORPTION_RATE,
	CAPTURE_RATE,
	FISSION_RATE,
	TRANSPORT_RATE,
	DIFFUSION_RATE,
	OUTSCATTER_RATE
} tallyType;


/* Tally spacing types */
typedef enum binSpacingTypes {
	EQUAL,
	LOGARITHMIC,
	OTHER
} binSpacingType;


/**
 * This class represents a set of tallies. A set of values
 * define the edges between bins for each tally. This class 
 * holds the edges, the centers between bins. It also allows 
 * for tallies to be made within each bin.
 */
class Tally{

protected:
	char* _tally_name;
	int _num_bins;
	double* _edges;
	double* _centers;
	double** _tallies;
	double _bin_delta;
	binSpacingType _bin_spacing;
	tallyDomainType _tally_domain;
	tallyType _tally_type;
    triggerType _trigger_type;
    float _trigger_precision;

	int _num_batches;
	double* _batch_mu;
	double* _batch_variance;
	double* _batch_std_dev;
	double* _batch_rel_err;
	bool _computed_statistics;

public:
	Tally(char* tally_name);
	virtual ~Tally();
	char* getTallyName();
	int getNumBins();
	double* getBinEdges();
	double* getBinCenters();
	double getBinDelta();
	double getBinDelta(double sample);
	binSpacingType getBinSpacingType();
	tallyDomainType getTallyDomainType();
	tallyType getTallyType();
	double** getTallies();
	double getTally(int bin_index, int batch_num);
	double getMaxTally();
	double getMinTally();
	int getBinIndex(double sample);

    double getMaxVariance();
    double getMaxStdDev();
    double getMaxRelErr();
	float getTriggerPrecision();
	triggerType getTriggerType();
    bool hasComputedBatchStatistics();

    /* IMPORTANT: The following five class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be transformed
     * into numpy arrays */
    void retrieveTallyEdges(double* data, int num_bins);
    void retrieveTallyCenters(double* data, int num_bins);
    void retrieveTallyMu(double* data, int num_bins);
    void retrieveTallyVariance(double* data, int num_bins);
    void retrieveTallyStdDev(double* data, int num_bins);
    void retrieveTallyRelErr(double* data, int num_bins);

	int getNumBatches();
	double* getBatchMu();
	double* getBatchVariance();
	double* getBatchStdDev();
	double* getBatchRelativeError();

    void setBinSpacingType(binSpacingType type);
	void setBinEdges(double* edges, int num_edges);
    void setPrecisionTrigger(triggerType trigger_type, float precision);
	void generateBinEdges(double start, double end, int num_bins, 
											binSpacingType type);
	void generateBinCenters();

	void setNumBatches(int num_batches);
    void incrementNumBatches(int num_batches);
    bool isPrecisionTriggered();

	void computeBatchStatistics();
	void computeScaledBatchStatistics(double scale_factor);
	void outputBatchStatistics(const char* filename);

	void tally(neutron* neutron, double weight);
	virtual void tally(neutron* neutron) =0;
//	virtual Tally* clone() =0;
};


class IsotopeTally: public Tally {

protected:
    Isotope* _isotope;
public:
	IsotopeTally(char* tally_name, Isotope* isotope)
			: Tally(tally_name){ 
				_tally_domain = ISOTOPE;
				_isotope = isotope; 
			}
	Isotope* getIsotope() { return _isotope; }
	virtual void tally(neutron* neutron) =0;
};


class MaterialTally: public Tally {

protected:
    Material* _material;
public:
	MaterialTally(char* tally_name, Material* material)
			: Tally(tally_name) { 
				_tally_domain = MATERIAL;
				_material = material; 
			}
	Material* getMaterial() { return _material; }
	virtual void tally(neutron* neutron) =0;
};


class RegionTally: public Tally {

protected:
    Region* _region;
public:
	RegionTally(char* tally_name, Region* region)
			: Tally(tally_name) { 
				_tally_domain = REGION;
				_region = region; 
			}
	Region* getRegion() { return _region; }
	virtual void tally(neutron* neutron) =0;
};


class GeometryTally: public Tally {

protected:
	Geometry* _geometry;
public:
	GeometryTally(char* tally_name, Geometry* geometry)
			: Tally(tally_name) {
				_tally_domain = GEOMETRY;
				_geometry = geometry;
			}
	Geometry* getGeometry() { return _geometry; }
	virtual void tally(neutron* neutron) =0;
};

/******************************************************************************/
/************************** Outscatter Rate Tallies ***************************/
/******************************************************************************/

class IsotopeOutscatterRateTally: public IsotopeTally {
public:
	IsotopeOutscatterRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope) {
				_tally_type = OUTSCATTER_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialOutscatterRateTally: public MaterialTally {

public:
	MaterialOutscatterRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = OUTSCATTER_RATE;
			}
	void tally(neutron* neutron);
};


class RegionOutscatterRateTally: public RegionTally {

public:
	RegionOutscatterRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = OUTSCATTER_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryOutscatterRateTally: public GeometryTally {

public:
	GeometryOutscatterRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = OUTSCATTER_RATE;
			}
	void tally(neutron* neutron);
};

/******************************************************************************/
/*************************** Collision Rate Tallies ***************************/
/******************************************************************************/

class IsotopeCollisionRateTally: public IsotopeTally {

public:
	IsotopeCollisionRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};


class MaterialCollisionRateTally: public MaterialTally {

public:
	MaterialCollisionRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};



class RegionCollisionRateTally: public RegionTally {

public:
	RegionCollisionRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};



class GeometryCollisionRateTally: public GeometryTally {

public:
	GeometryCollisionRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};



/******************************************************************************/
/**************************** Elastic Rate Tallies ****************************/
/******************************************************************************/
//FIXME - implement scattering matrix and override most Tally methods

class IsotopeElasticRateTally: public IsotopeTally {

public:
	IsotopeElasticRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialElasticRateTally: public MaterialTally {

public:
	MaterialElasticRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


class RegionElasticRateTally: public RegionTally {

public:
	RegionElasticRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryElasticRateTally: public GeometryTally {

public:
	GeometryElasticRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Absorption Rate Tallies ***************************/
/******************************************************************************/

class IsotopeAbsorptionRateTally: public IsotopeTally {
public:
	IsotopeAbsorptionRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialAbsorptionRateTally: public MaterialTally {

public:
	MaterialAbsorptionRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


class RegionAbsorptionRateTally: public RegionTally {

public:
	RegionAbsorptionRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryAbsorptionRateTally: public GeometryTally {

public:
	GeometryAbsorptionRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/**************************** Capture Rate Tallies ****************************/
/******************************************************************************/

class IsotopeCaptureRateTally: public IsotopeTally {

public:
	IsotopeCaptureRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialCaptureRateTally: public MaterialTally {

public:
	MaterialCaptureRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


class RegionCaptureRateTally: public RegionTally {

public:
	RegionCaptureRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryCaptureRateTally: public GeometryTally {

public:
	GeometryCaptureRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Fission Rate Tallies *****************************/
/******************************************************************************/

class IsotopeFissionRateTally: public IsotopeTally {

public:
	IsotopeFissionRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialFissionRateTally: public MaterialTally {

public:
	MaterialFissionRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


class RegionFissionRateTally: public RegionTally {

public:
	RegionFissionRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryFissionRateTally: public GeometryTally {

public:
	GeometryFissionRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Transport Rate Tallies ****************************/
/******************************************************************************/

class IsotopeTransportRateTally: public IsotopeTally {

public:
	IsotopeTransportRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialTransportRateTally: public MaterialTally {

public:
	MaterialTransportRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


class RegionTransportRateTally: public RegionTally {

public:
	RegionTransportRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryTransportRateTally: public GeometryTally {

public:
	GeometryTransportRateTally(char* tally_name, Geometry* geometry) 
			: GeometryTally(tally_name, geometry) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Diffusion Rate Tallies ****************************/
/******************************************************************************/

class IsotopeDiffusionRateTally: public IsotopeTally {

public:
	IsotopeDiffusionRateTally(char* tally_name, Isotope* isotope)
			: IsotopeTally(tally_name, isotope) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialDiffusionRateTally: public MaterialTally {

public:
	MaterialDiffusionRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


class RegionDiffusionRateTally: public RegionTally {

public:
	RegionDiffusionRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryDiffusionRateTally: public GeometryTally {

public:
	GeometryDiffusionRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Leakage Rate Tallies *****************************/
/******************************************************************************/
//FIXME

class MaterialLeakageRateTally: public MaterialTally {

public:
	MaterialLeakageRateTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = LEAKAGE_RATE;
			}
	void tally(neutron* neutron);
};


class RegionLeakageRateTally: public RegionTally {

public:
	RegionLeakageRateTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = LEAKAGE_RATE;
			}
	void tally(neutron* neutron);
};

class GeometryLeakageRateTally: public GeometryTally {

public:
	GeometryLeakageRateTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = LEAKAGE_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/****************************** Flux Tallies **********************************/
/******************************************************************************/

class MaterialFluxTally: public MaterialTally {

public:
	MaterialFluxTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material) {
				_tally_type = FLUX;
			}
	void tally(neutron* neutron);
};


class RegionFluxTally: public RegionTally {

public:
	RegionFluxTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region) {
				_tally_type = FLUX;
			}
	void tally(neutron* neutron);
};


class GeometryFluxTally: public GeometryTally {

public:
	GeometryFluxTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry) {
				_tally_type = FLUX;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Mean Lifetime Tallies ****************************/
/******************************************************************************/

class MaterialInterCollisionTimeTally: public MaterialTally {

public:
	MaterialInterCollisionTimeTally(char* tally_name, Material* material)
			: MaterialTally(tally_name, material){
 				_tally_type = INTERCOLLISION_TIME; 
			}
	void tally(neutron* neutron);
};



class RegionInterCollisionTimeTally: public RegionTally {

public:
	RegionInterCollisionTimeTally(char* tally_name, Region* region)
			: RegionTally(tally_name, region){
 				_tally_type = INTERCOLLISION_TIME; 
			}
	void tally(neutron* neutron);
};



class GeometryInterCollisionTimeTally: public GeometryTally {

public:
	GeometryInterCollisionTimeTally(char* tally_name, Geometry* geometry)
			: GeometryTally(tally_name, geometry){
 				_tally_type = INTERCOLLISION_TIME; 
			}
	void tally(neutron* neutron);
};


/**
 * Finds the bin index for a sample in a set of bins. If the samples
 * is outside the bounds of all bins, it returns infinity
 * @param sample the sample value of interest
 * @return the bin index for the sample
 */
inline int Tally::getBinIndex(double sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return a bin index for Tally %s since "
				 "the bins have not yet been created", _tally_name);

	/* Set index to infinity to begin with */
	int index = std::numeric_limits<int>::infinity();

	/* if the sample is equal to the last bin edge, return the last bin */
	if (sample == _edges[_num_bins])
		return _num_bins-1;

	/* Equally spaced bins */
	if (_bin_spacing == EQUAL)
		index = int((sample - _edges[0]) / _bin_delta);

	/* Logarithmically spaced bins */
	else if (_bin_spacing == LOGARITHMIC)
		index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

	/* If the bin_type == OTHER then the bin edges were not generated by
	 * generateEqualBinEdges, so use a binary search to find the bin */
	else
		index = findUpperIndex(_edges, _num_bins, 0, sample) - 1;

	/* If this sample was not contained within a bin set index to infinity*/
	if (index > _num_bins)
		index = std::numeric_limits<int>::infinity();

	return index;
}

#endif

#endif /* BINNER_H_ */
