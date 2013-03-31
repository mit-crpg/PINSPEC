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
#include <string.h>
#include <sstream>
#include "log.h"
#include "arraycreator.h"
#include "Neutron.h"
#include "Material.h"
#include "Isotope.h"
#include "Region.h"

class Geometry;
class DerivedTally;

#define NEUTRON_MASS 939565378           /* Mass of neutron in eV / c^2 */
#define LIGHT_SPEED 299792458               /* Speed of light in m / s */


/* The domain in which this tally resides */
typedef enum tallyDomainTypes {
	MATERIAL,
	ISOTOPE,
	REGION,
    GEOMETRY,
    UNDEFINED                       /* For derived type tallies */
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
/* DERIVED type tallies are returned when users apply mathematical operations 
 * to two or more tallies (ie, tally2 / tally2 or tally1 + tally2) */
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
    DERIVED
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
	Tally(const char* tally_name=(char*)"");
//    Tally(const Tally &tally);       /* Copy constructor - how to? */
//    Tally(const Tally *tally);             /* Copy constructor - how to? */
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

    double getMaxMu();
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

    void setTallyDomainType(tallyDomainType type);
    void setTallyType(tallyType type);
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
    void normalizeBatchMu();
	void outputBatchStatistics(const char* filename);
    void printTallies(bool uncertainties=false);
    Tally* clone();

	void tally(neutron* neutron, double weight);
	virtual void tally(neutron* neutron) =0;

    /* Perhaps we can generalize the following functions somewhat, but
     * this is what I came up with for now */
    DerivedTally* addIntegers(const int* amt, const int length);
    DerivedTally* addFloats(const float* amt, const int length);
    DerivedTally* addDoubles(const double* amt, const int length);

    DerivedTally* subtractIntegers(const int* amt, const int length);
    DerivedTally* subtractFloats(const float* amt, const int length);
    DerivedTally* subtractDoubles(const double* amt, const int length);

    DerivedTally* multiplyIntegers(const int* amt, const int length);
    DerivedTally* multiplyFloats(const float* amt, const int length);
    DerivedTally* multiplyDoubles(const double* amt, const int length);

    DerivedTally* divideIntegers(const int* amt, const int length);
    DerivedTally* divideFloats(const float* amt, const int length);
    DerivedTally* divideDoubles(const double* amt, const int length);

    /* Assignment operators - how to? */
//    Tally* operator=(Tally* tally);
//    Tally* operator=(const Tally& tally);

    /* Operator overloads */
    DerivedTally* operator+(Tally* tally);
    DerivedTally* operator-(Tally* tally);
    DerivedTally* operator*(Tally* tally);
    DerivedTally* operator/(Tally* tally);

    DerivedTally* operator+(const int amt);
    DerivedTally* operator-(const int amt);
    DerivedTally* operator*(const int amt);
    DerivedTally* operator/(const int amt);

    DerivedTally* operator+(const float amt);
    DerivedTally* operator-(const float amt);
    DerivedTally* operator*(const float amt);
    DerivedTally* operator/(const float amt);

    DerivedTally* operator+(const double amt);
    DerivedTally* operator-(const double amt);
    DerivedTally* operator*(const double amt);
    DerivedTally* operator/(const double amt);
};


class IsotopeTally: public Tally {

protected:
    Isotope* _isotope;
public:
	IsotopeTally(Isotope* isotope, const char* tally_name=(char*)"")
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
	MaterialTally(Material* material, const char* tally_name=(char*)"")
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
	RegionTally(Region* region, const char* tally_name=(char*)"")
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
	GeometryTally(Geometry* geometry, const char* tally_name=(char*)"")
			: Tally(tally_name) {
				_tally_domain = GEOMETRY;
				_geometry = geometry;
			}
	Geometry* getGeometry() { return _geometry; }
	virtual void tally(neutron* neutron) =0;
};


/******************************************************************************/
/*************************** Collision Rate Tallies ***************************/
/******************************************************************************/

class IsotopeCollisionRateTally: public IsotopeTally {

public:
	IsotopeCollisionRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};


class MaterialCollisionRateTally: public MaterialTally {

public:
	MaterialCollisionRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};



class RegionCollisionRateTally: public RegionTally {

public:
	RegionCollisionRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name){
 				_tally_type = COLLISION_RATE; 
			}
	void tally(neutron* neutron);
};



class GeometryCollisionRateTally: public GeometryTally {

public:
	GeometryCollisionRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name){
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
	IsotopeElasticRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialElasticRateTally: public MaterialTally {

public:
	MaterialElasticRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


class RegionElasticRateTally: public RegionTally {

public:
	RegionElasticRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryElasticRateTally: public GeometryTally {

public:
	GeometryElasticRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name){
				_tally_type = ELASTIC_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Absorption Rate Tallies ***************************/
/******************************************************************************/

class IsotopeAbsorptionRateTally: public IsotopeTally {

public:
	IsotopeAbsorptionRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialAbsorptionRateTally: public MaterialTally {

public:
	MaterialAbsorptionRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


class RegionAbsorptionRateTally: public RegionTally {

public:
	RegionAbsorptionRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryAbsorptionRateTally: public GeometryTally {

public:
	GeometryAbsorptionRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name) {
				_tally_type = ABSORPTION_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/**************************** Capture Rate Tallies ****************************/
/******************************************************************************/

class IsotopeCaptureRateTally: public IsotopeTally {

public:
	IsotopeCaptureRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialCaptureRateTally: public MaterialTally {

public:
	MaterialCaptureRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


class RegionCaptureRateTally: public RegionTally {

public:
	RegionCaptureRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryCaptureRateTally: public GeometryTally {

public:
	GeometryCaptureRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name) {
				_tally_type = CAPTURE_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Fission Rate Tallies *****************************/
/******************************************************************************/

class IsotopeFissionRateTally: public IsotopeTally {

public:
	IsotopeFissionRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialFissionRateTally: public MaterialTally {

public:
	MaterialFissionRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


class RegionFissionRateTally: public RegionTally {

public:
	RegionFissionRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryFissionRateTally: public GeometryTally {

public:
	GeometryFissionRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name) {
				_tally_type = FISSION_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Transport Rate Tallies ****************************/
/******************************************************************************/

class IsotopeTransportRateTally: public IsotopeTally {

public:
	IsotopeTransportRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialTransportRateTally: public MaterialTally {

public:
	MaterialTransportRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


class RegionTransportRateTally: public RegionTally {

public:
	RegionTransportRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryTransportRateTally: public GeometryTally {

public:
	GeometryTransportRateTally(Geometry* geometry, const char* tally_name=(char*)"") 
			: GeometryTally(geometry, tally_name) {
				_tally_type = TRANSPORT_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Diffusion Rate Tallies ****************************/
/******************************************************************************/

class IsotopeDiffusionRateTally: public IsotopeTally {

public:
	IsotopeDiffusionRateTally(Isotope* isotope, const char* tally_name=(char*)"")
			: IsotopeTally(isotope, tally_name) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


class MaterialDiffusionRateTally: public MaterialTally {

public:
	MaterialDiffusionRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


class RegionDiffusionRateTally: public RegionTally {

public:
	RegionDiffusionRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


class GeometryDiffusionRateTally: public GeometryTally {

public:
	GeometryDiffusionRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name) {
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
	MaterialLeakageRateTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = LEAKAGE_RATE;
			}
	void tally(neutron* neutron);
};


class RegionLeakageRateTally: public RegionTally {

public:
	RegionLeakageRateTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = LEAKAGE_RATE;
			}
	void tally(neutron* neutron);
};

class GeometryLeakageRateTally: public GeometryTally {

public:
	GeometryLeakageRateTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name) {
				_tally_type = LEAKAGE_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/****************************** Flux Tallies **********************************/
/******************************************************************************/

class MaterialFluxTally: public MaterialTally {

public:
	MaterialFluxTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name) {
				_tally_type = FLUX;
			}
	void tally(neutron* neutron);
};


class RegionFluxTally: public RegionTally {

public:
	RegionFluxTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name) {
				_tally_type = FLUX;
			}
	void tally(neutron* neutron);
};


class GeometryFluxTally: public GeometryTally {

public:
	GeometryFluxTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name) {
				_tally_type = FLUX;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Mean Lifetime Tallies ****************************/
/******************************************************************************/

class MaterialInterCollisionTimeTally: public MaterialTally {

public:
	MaterialInterCollisionTimeTally(Material* material, const char* tally_name=(char*)"")
			: MaterialTally(material, tally_name){
 				_tally_type = INTERCOLLISION_TIME; 
			}
	void tally(neutron* neutron);
};



class RegionInterCollisionTimeTally: public RegionTally {

public:
	RegionInterCollisionTimeTally(Region* region, const char* tally_name=(char*)"")
			: RegionTally(region, tally_name){
 				_tally_type = INTERCOLLISION_TIME; 
			}
	void tally(neutron* neutron);
};



class GeometryInterCollisionTimeTally: public GeometryTally {

public:
	GeometryInterCollisionTimeTally(Geometry* geometry, const char* tally_name=(char*)"")
			: GeometryTally(geometry, tally_name){
 				_tally_type = INTERCOLLISION_TIME; 
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/****************************** Derived Tallies *******************************/
/******************************************************************************/

class DerivedTally: public Tally {

public:
	DerivedTally(const char* tally_name=(char*)"")
			: Tally(tally_name){ 
				_tally_domain = UNDEFINED;
                _tally_type = DERIVED;
			}
	void tally(neutron* neutron);
    void setTallyName(char* tally_name);
    void setTallies(double** tallies);
    void setBatchMu(double* batch_mu);
    void setBatchVariance(double* batch_variance);
    void setBatchStdDev(double* batch_std_dev);
    void setBatchRelErr(double* batch_rel_err);
    void setComputedBatchStatistics(bool computed);
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

	/* Logarithmically spaced bins */
	if (_bin_spacing == LOGARITHMIC)
		index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

	/* Equally spaced bins */
	else if (_bin_spacing == EQUAL)
		index = int((sample - _edges[0]) / _bin_delta);

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
