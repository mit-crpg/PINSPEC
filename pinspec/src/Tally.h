/**
 * @file Tally.h
 * @brief The Tally superclass and subclasses.
 * @author William Boyd (wboyd.mit.edu)
 * @date March 4, 2013
 */

#ifndef TALLY_H_
#define TALLY_H_

#ifdef __cplusplus
#include <limits>
#include <cmath>
#include <math.h>
#include <string.h>
#include <sstream>
#include "log.h"
#include "arraycreator.h"
#include "Neutron.h"
#include "Region.h"
#include "Material.h"
#include "Isotope.h"
#include "Region.h"
#endif

class Geometry;
class DerivedTally;

/** Mass of a neutron in \f$ \frac{eV}{c^2} \f$ */
#define NEUTRON_MASS 939565378

/** The speed of light in meters per second */
#define LIGHT_SPEED 299792458


/**
 * @enum tallyDomainTypes
 * @brief The types of domains in which tally may reside.
 */

/**
 * @var tallyDomainType
 * @brief A domain type within which a tally may reside.
 */
typedef enum tallyDomainTypes {
    /** A tally may reside within a material */
    MATERIAL,
    /** A tally may reside within an isotope */
    ISOTOPE,
    /** A tally may reside within a region */
    REGION,
    /** A tally may reside within the entire geometry */
    GEOMETRY,
    /** A tally may have an undefined domain - only used for DERIVED
     * tally types */
    UNDEFINED
} tallyDomainType;


/**
 * @enum triggerTypes
 * @brief The type of precision trigger for a tally object.
 * @details Precision triggers represent a threshold on the level of precision
 *          the user requires before ending the simulation.
 */

/**
 * @var triggerType
 * @brief The type of precision trigger for a tally object.
 * @details Precision triggers represent a threshold on the level of precision
 *          the user requires before ending the simulation.
 */
typedef enum triggerTypes {
    /** A precision trigger on the maximum variance for a tally */
    VARIANCE,
    /** A precision trigger on the maximum standard deviation for a tally */
    STANDARD_DEVIATION,
    /** A precision trigger on the maximum relative error for a tally */
    RELATIVE_ERROR,
    /** The absence of a precision trigger for a tally */
    NONE
} triggerType;


/* Type of tallies */
/* DERIVED type tallies are returned when users apply mathematical operations 
 * to two or more tallies (ie, tally2 / tally2 or tally1 + tally2) */
/**
 * @enum tallyTypes
 * @brief The types of tallies which may be instantiated and used in a 
 *        PINSPEC simulation.
 */

/**
 * @var tallyType
 * @brief The types of tallies which may be instantiated and used in a 
 *        PINSPEC simulation.
 */
typedef enum tallyTypes {
    /** A tally of the flux */
    FLUX,
    /** A tally of the leakage rate */
    LEAKAGE_RATE,
    /** A tally of the collision rate */
    COLLISION_RATE,
    /** A tally of the intercollision time */
    INTERCOLLISION_TIME,
    /** A tally of the elastic scattering rate */
    ELASTIC_RATE,
    /** A tally of the group-to-group scattering rate */
    GROUP_TO_GROUP_RATE,
    /** A tally of the out-scattering reaction rate */
    OUTSCATTER_RATE,
    /** A tally of the absorption rate */
    ABSORPTION_RATE,
    /** A tally of the capture rate */
    CAPTURE_RATE,
    /** A tally of the fission rate */
    FISSION_RATE,
    /** A tally of the transport rate */
    TRANSPORT_RATE,
    /** A tally of the diffusion rate */
    DIFFUSION_RATE,
    /** A derived tally type - used to define the tally returned from tally 
     * arithmetic such as would be the case in Python:
     *
     * @code
     *         derived_tally = tally1 + tally2
     * @endcode
     *
     */
    DERIVED
} tallyType;



/** 
 * @enum binSpacingTypes
 * @brief The spacing between bin edges for tallies
 */

/**
 * @var binSpacingType
 * @brief The spacing between bin edges for a tally
 */
typedef enum binSpacingTypes {
    /** Equally spaced bins */
    EQUAL,
    /** Logarithmically spaced bins */
    LOGARITHMIC,
    /** Bin edges without a specified pattern */
    OTHER
} binSpacingType;


/**
 * This class represents a set of tallies. A set of values
 * define the edges between bins for each tally. This class 
 * holds the edges, the centers between bins. It also allows 
 * for tallies to be made within each bin.
 */
/**
 * @class Tally Tally.h "pinspec/src/Tally.h"
 * @brief A Tally reprsents a set of bins for tallying some quantity.
 * @details A set of values define the edges between bins for each tally. This
 *          class holds the edges, the centers between bins and the tallies
 *          within each bin. The Tally class knows how to compute batch-based
 *          statistics for each tally bin. The Tally class is an abstract class
 *          and must be implemented for each specific type of Tally the developer
 *          might wish to define.
 */
class Tally {

protected:
    /** The user-specified name for the tally */
    char* _tally_name;
    /** The number of tally bins */
    int _num_bins;
    /** The number of tally edges */
    int _num_edges;
    /** The array of bin edges between tally bins */
    double* _edges;
    /** The array of bin center values */
    double* _centers;
    /** A 2D array of tallies for each bin and each batch */
    double** _tallies;
    /** Equal / logarithmic spacing between bins if defined on a uniform grid */
    double _bin_delta;
    /** The spacing type between bins */
    binSpacingType _bin_spacing;
    /** The domain in which this tally resides */
    tallyDomainType _tally_domain;
    /** The type of tally */
    tallyType _tally_type;
    /** The precision trigger for this tally */
    triggerType _trigger_type;
    /** The trigger precision for this tally */
    float _trigger_precision;

    /** The number of batches in the PINSPEC simulation */
    int _num_batches;
    /** The batch average for each tally bin */
    double* _batch_mu;
    /** The batch variance for each tally bin */
    double* _batch_variance;
    /** The batch standard deviation for each tally bin */
    double* _batch_std_dev;
    /** The batch relative error for each tally bin */
    double* _batch_rel_err;
    /** Whether or not batch statistics have been computed */
    bool _computed_statistics;
    /** Whether or not bin size has be squared for group-to-group xs */
    bool _group_expand_bins;

public:
    Tally(char* tally_name=(char*)"");
    virtual ~Tally();
    char* getTallyName();
    int getNumBins();
    int getNumEdges();
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
    bool hasExpandedGroupBins();

    /* IMPORTANT: The following six class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be 
     * transformed into numpy arrays */
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
    void setGroupExpandBins(bool expand_bins);
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
    void tallyGroup(neutron* neutron, double weight);

    /**
     * @brief A virtual method to tally a neutron which must be implemented by
     *        Tally subclasses.
     * @param neutron the neutron of interest
     */
    virtual void tally(neutron* neutron) =0;

    DerivedTally* tile(const int num_tiles);

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

    /* Operator overloads for tally arithmetic */
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


/**
 * @class IsotopeTally Tally.h "pinspec/src/Tally.h"
 * @brief An abstract class for tallies with ISOTOPE domain type.
 * @details An IsotopeTally is for tallying a specific reaction rate for
 *          an isotope within a material, region or the geometry. Isotope
 *          tallies are especially useful for computing resonance integrals.
 *          The IsotopeTally is an abstract class and must be implemented for
 *          each tallyType.
 */


class IsotopeTally: public Tally {

protected:
    /** A pointer to the isotope for this tally */
    Isotope* _isotope;

public:
    /** 
     * @brief IsotopeTally constructor calls the Tally constructor and sets
     *        the tally domain to ISOTOPE.
     * @param isotope a pointer to the isotope to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeTally(Isotope* isotope, char* tally_name=(char*)"")
      : Tally(tally_name){ 
          _tally_domain = ISOTOPE;
	  _isotope = isotope; 
    }

    
    /**
     * @brief Returns the isotope for this tally.
     * @return a pointer to the isotope
     */
    Isotope* getIsotope() { return _isotope; }

    /** @brief A method to tally a neutron to be implemented by subclasses. */
    virtual void tally(neutron* neutron) =0;
};


/**
 * @class MaterialTally Tally.h "pinspec/src/Tally.h"
 * @brief An abstract class for tallies with MATERIAL domain type.
 * @details A MaterialTally is for tallying within a material. The 
 *          MaterialTally is an abstract class and must be implemented for
 *          each tallyType.
 */
class MaterialTally: public Tally {

protected:
    /** A pointer to the material in which this tally resides */
    Material* _material;
public:
    /** 
     * @brief MaterialTally constructor calls the Tally constructor and sets
     *        the tally domain to MATERIAL.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialTally(Material* material, char* tally_name=(char*)"")
      : Tally(tally_name) { 
          _tally_domain = MATERIAL;
	  _material = material; 
    }

    /**
     * @brief Returns the material in which this tally resides.
     * @return a pointer to the material
     */
    Material* getMaterial() { return _material; }
     
    /** @brief A method to tally a neutron to be implemented by subclasses. */
    virtual void tally(neutron* neutron) =0;
};


/**
 * @class RegionTally Tally.h "pinspec/src/Tally.h"
 * @brief An abstract class for tallies with REGION domain type.
 * @details A RegionTally is for tallying within a region. The 
 *          RegionTally is an abstract class and must be implemented for
 *          each tallyType.
 */
class RegionTally: public Tally {

protected:
    /** A pointer to the region in which this tally resides */
    Region* _region;
public:
    /** 
     * @brief RegionTally constructor calls the Tally constructor and sets
     *        the tally domain to REGION.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionTally(Region* region, char* tally_name=(char*)"")
      : Tally(tally_name) { 
          _tally_domain = REGION;
	  _region = region; 
    }

    /**
     * @brief Returns the region in which this tally resides.
     * @return a pointer to the region
     */
    Region* getRegion() { return _region; }

    /** @brief A method to tally a neutron to be implemented by subclasses. */
    virtual void tally(neutron* neutron) =0;
};


/**
 * @class GeometryTally Tally.h "pinspec/src/Tally.h"
 * @brief An abstract class for tallies with MATERIAL domain type.
 * @details A GeometryTally is for tallying within a geometry. The 
 *          GeometryTally is an abstract class and must be implemented for
 *          each tallyType.
 */
class GeometryTally: public Tally {

protected:
    /** A pointer to the geometry in which this tally resides */
    Geometry* _geometry;
public:
    /** 
     * @brief GeometryTally constructor calls the Tally constructor and sets
     *        the tally domain to GEOMETRY.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryTally(Geometry* geometry, char* tally_name=(char*)"")
      : Tally(tally_name) {
          _tally_domain = GEOMETRY;
	  _geometry = geometry;
    }

    /**
     * @brief Returns the geometry in which this tally resides.
     * @return a pointer to the geometry
     */
    Geometry* getGeometry() { return _geometry; }
     
    /** @brief A method to tally a neutron to be implemented by subclasses. */
    virtual void tally(neutron* neutron) =0;
};


/******************************************************************************/
/*************************** Collision Rate Tallies ***************************/
/******************************************************************************/

/**
 * @class IsotopeCollisionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the collision rate for an isotope.
 */
class IsotopeCollisionRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeCollisionRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to COLLISION_RATE.
     * @param isotope a pointer to the isotope within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeCollisionRateTally(Isotope* isotope, char* tally_name=(char*)"")
        : IsotopeTally(isotope, tally_name){
            _tally_type = COLLISION_RATE; 
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialCollisionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the collision rate within a material.
 */
class MaterialCollisionRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialCollisionRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to COLLISION_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialCollisionRateTally(Material* material, 
			       char* tally_name=(char*)"")
        : MaterialTally(material, tally_name){
            _tally_type = COLLISION_RATE; 
    }

    void tally(neutron* neutron);
};



/**
 * @class RegionCollisionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the collision rate within a region.
 */
class RegionCollisionRateTally: public RegionTally {

public:
    /** 
     * @brief RegionCollisionRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to COLLISION_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionCollisionRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name){
            _tally_type = COLLISION_RATE; 
    }

    void tally(neutron* neutron);
};



/**
 * @class GeometryCollisionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the collision rate within the geometry.
 */
class GeometryCollisionRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryCollisionRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to COLLISION_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryCollisionRateTally(Geometry* geometry, 
			       char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name){
             _tally_type = COLLISION_RATE; 
    }

    void tally(neutron* neutron);
};



/******************************************************************************/
/**************************** Elastic Rate Tallies ****************************/
/******************************************************************************/
//FIXME - implement scattering matrix and override most Tally methods

/**
 * @class IsotopeElasticRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the elastic scattering rate for an isotope.
 */
class IsotopeElasticRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeElasticRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to ELASTIC_RATE.
     * @param isotope a pointer to the isotope for which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeElasticRateTally(Isotope* isotope, char* tally_name=(char*)"")
        : IsotopeTally(isotope, tally_name){
            _tally_type = ELASTIC_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialElasticRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the elastic scattering rate within a material.
 */
class MaterialElasticRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialCollisionRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to ELASTIC_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialElasticRateTally(Material* material, 
			     char* tally_name=(char*)"")
        : MaterialTally(material, tally_name){
            _tally_type = ELASTIC_RATE;
    }
  
    void tally(neutron* neutron);
};


/**
 * @class RegionElasticRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the elastic scattering rate within a region.
 */
class RegionElasticRateTally: public RegionTally {

public:
    /** 
     * @brief RegionElasticRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to ELASTIC_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionElasticRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name){
            _tally_type = ELASTIC_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryElasticRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the elastic scattering rate within the geometry.
 */
class GeometryElasticRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryElasticRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to ELASTIC_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryElasticRateTally(Geometry* geometry, 
			     char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name){
            _tally_type = ELASTIC_RATE;
    }

    void tally(neutron* neutron);
};

/******************************************************************************/
/**************************** GROUP Rate Tallies ****************************/
/******************************************************************************/
//FIXME - implement scattering matrix and override most Tally methods

/**
 * @class IsotopeGroupRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the group-to-group scattering rate for an isotope.
 */
class IsotopeGroupRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeGroupRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to GROUP_TO_GROUP_RATE.
     * @param isotope a pointer to the isotope for which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeGroupRateTally(Isotope* isotope, char* tally_name=(char*)"")
        : IsotopeTally(isotope, tally_name){
            _tally_type = GROUP_TO_GROUP_RATE;
            _group_expand_bins = false;
    }

    void tally(neutron* neutron);
};

/**
 * @class MaterialGroupRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the group-to-group scattering rate within a material.
 */
class MaterialGroupRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialGroupRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to GROUP_TO_GROUP_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialGroupRateTally(Material* material, 
			   char* tally_name=(char*)"")
        : MaterialTally(material, tally_name){
            _tally_type = GROUP_TO_GROUP_RATE;
            _group_expand_bins = false;
    }
  
    void tally(neutron* neutron);
};


/**
 * @class RegionGroupRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the group-to-group scattering rate within a region.
 */
class RegionGroupRateTally: public RegionTally {

public:
    /** 
     * @brief RegionGroupRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to GROUP_TO_GROUP_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionGroupRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name){
            _tally_type = GROUP_TO_GROUP_RATE;
            _group_expand_bins = false;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryGroupRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the group-to-group scattering rate within the geometry.
 */
class GeometryGroupRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryGroupRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to GROUP_TO_GROUP_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryGroupRateTally(Geometry* geometry, 
			   char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name){
            _tally_type = GROUP_TO_GROUP_RATE;
            _group_expand_bins = false;
    }

    void tally(neutron* neutron);
};

/******************************************************************************/
/************************ Outter Scatter Rate Tallies  ************************/
/******************************************************************************/
/**
 * @class IsotopeOutScatterRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the outter scattering rate for an isotope.
 */
class IsotopeOutScatterRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeOutScatterRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to OUTSCATTER_RATE.
     * @param isotope a pointer to the isotope for which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeOutScatterRateTally(Isotope* isotope, char* tally_name=(char*)"")
        : IsotopeTally(isotope, tally_name){
            _tally_type = OUTSCATTER_RATE;
    }

    void tally(neutron* neutron);
};

/**
 * @class MaterialOutScatterRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the out scattering rate within a material.
 */
class MaterialOutScatterRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialOutScatterRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to OUTSCATTER_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialOutScatterRateTally(Material* material, 
			   char* tally_name=(char*)"")
        : MaterialTally(material, tally_name){
            _tally_type = OUTSCATTER_RATE;
    }
  
    void tally(neutron* neutron);
};


/**
 * @class RegionOutScatterRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the outter scattering rate within a region.
 */
class RegionOutScatterRateTally: public RegionTally {

public:
    /** 
     * @brief RegionOutScatterRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to OUTSCATTER_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionOutScatterRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name){
	_tally_type = OUTSCATTER_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryOutScatterRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the outter scattering rate within the geometry.
 */
class GeometryOutScatterRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryOutScatterRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to OUTSCATTER_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryOutScatterRateTally(Geometry* geometry, 
			   char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name){
	_tally_type = OUTSCATTER_RATE;
    }

    void tally(neutron* neutron);
};

/******************************************************************************/
/************************** Absorption Rate Tallies ***************************/
/******************************************************************************/

/**
 * @class IsotopeAbsorptionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the absorption rate for an isotope.
 */
class IsotopeAbsorptionRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeAbsorptionRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to ABSORPTION_RATE.
     * @param isotope  a pointer to the isotope within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeAbsorptionRateTally(Isotope* isotope, 
			       char* tally_name=(char*)"")
        : IsotopeTally(isotope, tally_name) {
            _tally_type = ABSORPTION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialAbsorptionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the absorption rate within a material
 */
class MaterialAbsorptionRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialAbsorptionRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to ABSORPTION_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialAbsorptionRateTally(Material* material, 
				char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
            _tally_type = ABSORPTION_RATE;
    }
  
    void tally(neutron* neutron);
};


/**
 * @class RegionAbsorptionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the absorption rate within a region.
 */
class RegionAbsorptionRateTally: public RegionTally {

public:
    /** 
     * @brief RegionAbsorptionRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to ABSORPTION_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionAbsorptionRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
            _tally_type = ABSORPTION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryAbsorptionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the absorption within the geometry.
 */
class GeometryAbsorptionRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryAbsorptionRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to ABSORPTION_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryAbsorptionRateTally(Geometry* geometry, 
				char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name) {
            _tally_type = ABSORPTION_RATE;
    }

    void tally(neutron* neutron);
};


/******************************************************************************/
/**************************** Capture Rate Tallies ****************************/
/******************************************************************************/


/**
 * @class IsotopeCaptureRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the capture rate for an isotope.
 */
class IsotopeCaptureRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeCaptureRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to CAPTURE_RATE.
     * @param isotope a pointer to the isotope within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeCaptureRateTally(Isotope* isotope, char* tally_name=(char*)"")
	: IsotopeTally(isotope, tally_name) {
            _tally_type = CAPTURE_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialCaptureRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the collision rate within a material.
 */
class MaterialCaptureRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialCaptureRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to CAPTURE_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialCaptureRateTally(Material* material, 
			     char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
            _tally_type = CAPTURE_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class RegionCaptureRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the capture rate within a region.
 */
class RegionCaptureRateTally: public RegionTally {

public:
    /** 
     * @brief RegionCaptureRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to CAPTURE_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionCaptureRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
            _tally_type = CAPTURE_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryCaptureRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the capture rate within the geometry.
 */
class GeometryCaptureRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryCaptureRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to CAPTURE_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryCaptureRateTally(Geometry* geometry, 
			     char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name) {
            _tally_type = CAPTURE_RATE;
    }

    void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Fission Rate Tallies *****************************/
/******************************************************************************/


/**
 * @class IsotopeFissionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the fission rate for an isotope.
 */
class IsotopeFissionRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeFissionRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to FISSION_RATE.
     * @param isotope a pointer to the isotope within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeFissionRateTally(Isotope* isotope, char* tally_name=(char*)"")
	: IsotopeTally(isotope, tally_name) {
            _tally_type = FISSION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialFissionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the capture rate within a material.
 */
class MaterialFissionRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialFissionRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to FISSION_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialFissionRateTally(Material* material, 
			     char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
            _tally_type = FISSION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class RegionFissionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the fission rate within a region..
 */
class RegionFissionRateTally: public RegionTally {

public:
    /** 
     * @brief RegionFissionRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to FISSION_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionFissionRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
            _tally_type = FISSION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryFissionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the capture rate for an isotope.
 */
class GeometryFissionRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryFissionRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to FISSION_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryFissionRateTally(Geometry* geometry, 
			     char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name) {
	    _tally_type = FISSION_RATE;
    }

    void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Transport Rate Tallies ****************************/
/******************************************************************************/

/**
 * @class IsotopeTransportRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the transport rate for an isotope.
 */
class IsotopeTransportRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeTransportRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to TRANSPORT_RATE.
     * @param isotope a pointer to the isotope within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeTransportRateTally(Isotope* isotope, char* tally_name=(char*)"")
        : IsotopeTally(isotope, tally_name) {
            _tally_type = TRANSPORT_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialTransportRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the transport rate within a material.
 */
class MaterialTransportRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialTransportRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to TRANSPORT_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialTransportRateTally(Material* material, 
			       char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
            _tally_type = TRANSPORT_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class RegionTransportRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the transport rate within a region.
 */
class RegionTransportRateTally: public RegionTally {

public:
    /** 
     * @brief RegionTransportRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to TRANSPORT_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionTransportRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
	    _tally_type = TRANSPORT_RATE;
    }
	
    void tally(neutron* neutron);
};


/**
 * @class GeometryTransportRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the transport rate within the geometry.
 */
class GeometryTransportRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryTransportRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to TRANSPORT_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryTransportRateTally(Geometry* geometry, 
			       char* tally_name=(char*)"") 
        : GeometryTally(geometry, tally_name) {
            _tally_type = TRANSPORT_RATE;
    }

    void tally(neutron* neutron);
};


/******************************************************************************/
/************************** Diffusion Rate Tallies ****************************/
/******************************************************************************/

/**
 * @class IsotopeDiffusionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the diffusion rate for an isotope.
 */
class IsotopeDiffusionRateTally: public IsotopeTally {

public:
    /** 
     * @brief IsotopeDiffusionRateTally constructor calls the IsotopeTally
     *        constructor and Tally constructors and sets the tally type
     *        to DIFFUSION_RATE.
     * @param isotope a pointer to the isotope within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    IsotopeDiffusionRateTally(Isotope* isotope, char* tally_name=(char*)"")
	: IsotopeTally(isotope, tally_name) {
            _tally_type = DIFFUSION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class MaterialDiffusionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the diffusion rate within a material.
 */
class MaterialDiffusionRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialDiffusionRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to DIFFUSION_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialDiffusionRateTally(Material* material, 
			       char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
            _tally_type = DIFFUSION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class RegionDiffusionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the diffusion rate within a region.
 */
class RegionDiffusionRateTally: public RegionTally {

public:
    /** 
     * @brief RegionDiffusionRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to DIFFUSION_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionDiffusionRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
            _tally_type = DIFFUSION_RATE;
    }

    void tally(neutron* neutron);
};


/**
 * @class GeometryDiffusionRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the capture rate for an isotope.
 */
class GeometryDiffusionRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryDiffusionRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to DIFFUSION_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryDiffusionRateTally(Geometry* geometry, 
			       char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name) {
            _tally_type = DIFFUSION_RATE;
			}
	void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Leakage Rate Tallies *****************************/
/******************************************************************************/
//FIXME

/**
 * @class MaterialLeakageRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the leakage rate for the regions filled by
 *        a material.
 */
class MaterialLeakageRateTally: public MaterialTally {

public:
    /** 
     * @brief MaterialLeakageRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to LEAKAGE_RATE.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialLeakageRateTally(Material* material, 
			     char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
            _tally_type = LEAKAGE_RATE;
    }
  
    void tally(neutron* neutron);
};


/**
 * @class RegionLeakageRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the leakage rate for a region.
 */
class RegionLeakageRateTally: public RegionTally {

public:
    /** 
     * @brief RegionLeakageRateTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to LEAKAGE_RATE.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionLeakageRateTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
             _tally_type = LEAKAGE_RATE;
    }
  
    void tally(neutron* neutron);
};


/**
 * @class GeometryLeakageRateTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the leakage rate for the geometry.
 */
class GeometryLeakageRateTally: public GeometryTally {

public:
    /** 
     * @brief GeometryLeakageRateTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to LEAKAGE_RATE.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
     GeometryLeakageRateTally(Geometry* geometry, 
			      char* tally_name=(char*)"")
         : GeometryTally(geometry, tally_name) {
	     _tally_type = LEAKAGE_RATE;
    }
	
    void tally(neutron* neutron);
};


/******************************************************************************/
/****************************** Flux Tallies **********************************/
/******************************************************************************/

/**
 * @class MaterialFluxTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the flux in all regions filled by a material.
 */
class MaterialFluxTally: public MaterialTally {

public:
    /** 
     * @brief IsotopeCaptureRateTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to FLUX.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialFluxTally(Material* material, char* tally_name=(char*)"")
        : MaterialTally(material, tally_name) {
	    _tally_type = FLUX;
    }
       
    void tally(neutron* neutron);
};


/**
 * @class RegionFluxTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the flux for a region.
 */
class RegionFluxTally: public RegionTally {

public:
    /** 
     * @brief RegionFluxTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to FLUX.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionFluxTally(Region* region, char* tally_name=(char*)"")
        : RegionTally(region, tally_name) {
            _tally_type = FLUX;
      }

    void tally(neutron* neutron);
};


/**
 * @class GeometryFluxTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the flux for the geometry.
 */
class GeometryFluxTally: public GeometryTally {

public:
    /** 
     * @brief GeometryFluxTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to FLUX.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryFluxTally(Geometry* geometry, char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name) {
            _tally_type = FLUX;
    }

    void tally(neutron* neutron);
};


/******************************************************************************/
/*************************** Mean Lifetime Tallies ****************************/
/******************************************************************************/


/**
 * @class MaterialInterCollisionTimeTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the time between collision with a material.
 */
class MaterialInterCollisionTimeTally: public MaterialTally {

public:
    /** 
     * @brief MaterialInterCollisionTimeTally constructor calls the MaterialTally
     *        constructor and Tally constructors and sets the tally type
     *        to INTERCOLLISION_TIME.
     * @param material a pointer to the material within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    MaterialInterCollisionTimeTally(Material* material, 
				    char* tally_name=(char*)"")
        : MaterialTally(material, tally_name){
            _tally_type = INTERCOLLISION_TIME; 
    }
  
    void tally(neutron* neutron);
};


/**
 * @class RegionInterCollisionTimeTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the time between collisions for the region.
 */
class RegionInterCollisionTimeTally: public RegionTally {

public:
    /** 
     * @brief RegionInterCollisionTimeTally constructor calls the RegionTally
     *        constructor and Tally constructors and sets the tally type
     *        to INTERCOLLISION_TIME.
     * @param region a pointer to the region within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    RegionInterCollisionTimeTally(Region* region, 
				  char* tally_name=(char*)"")
        : RegionTally(region, tally_name){
            _tally_type = INTERCOLLISION_TIME; 
    }
  
    void tally(neutron* neutron);
};


/**
 * @class GeometryInterCollisionTimeTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallying the time between collisions for the geometry.
 */
class GeometryInterCollisionTimeTally: public GeometryTally {

public:
    /** 
     * @brief GeometryInterCollisionTimeTally constructor calls the GeometryTally
     *        constructor and Tally constructors and sets the tally type
     *        to INTERCOLLISION_TIME.
     * @param geometry a pointer to the geometry within which to tally
     * @param tally_name a character array for the tally name (optional)
     */
    GeometryInterCollisionTimeTally(Geometry* geometry, 
				    char* tally_name=(char*)"")
        : GeometryTally(geometry, tally_name){
            _tally_type = INTERCOLLISION_TIME; 
    }

    void tally(neutron* neutron);
};


/******************************************************************************/
/****************************** Derived Tallies *******************************/
/******************************************************************************/

/**
 * @class DerivedTally Tally.h "pinspec/src/Tally.h"
 * @brief A class for tallies resulting from tally arithmetic operations.
 */
class DerivedTally: public Tally {

public: 
    /** 
     * @brief DerivedTally constructor calls the Tally constructor and sets 
     *        the tally domain type to UNDEFINED and the tally type to 
     *        DERIVED.
     * @param tally_name a character array for the tally name (optional)
     */
    DerivedTally(char* tally_name=(char*)"")
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
 * @brief Finds the bin index for a sample in a set of bins. If the samples
 *        is outside the bounds of all bins, it returns infinity
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
    if (sample == _edges[_num_edges - 1])
        return _num_edges - 2;

    /* Logarithmically spaced bins */
    if (_bin_spacing == LOGARITHMIC)
        index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

    /* Equally spaced bins */
    else if (_bin_spacing == EQUAL)
        index = int((sample - _edges[0]) / _bin_delta);

    /* If the bin_type == OTHER then the bin edges were not generated by
     * generateEqualBinEdges, so use a binary search to find the bin */
    else
        index = findUpperIndex(_edges, _num_edges - 1, 0, sample) - 1;

    /* If this sample was not contained within a bin set index to infinity*/
    if (index > _num_edges - 1)
        index = std::numeric_limits<int>::infinity();

    return index;
}


#endif /* TALLY_H_ */
