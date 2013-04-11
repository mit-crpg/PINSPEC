/**
 * @file TallyBank.h
 * @brief The TallyBank static class.
 * @author William Boyd (wboyd.mit.edu)
 * @date March 20, 2013
 */

#ifndef TALLYBANK_H_
#define TALLYBANK_H_

#ifdef __cplusplus
#include <set>
#include <map>
#include <utility>
#include <string.h>
#include <sstream>
#include "Tally.h"
#include "Region.h"
#include "Material.h"
#include "Isotope.h"
#include "Geometry.h"
#endif

#ifndef TALLYBANK_C
    /** The filename for this tally's output file of batch-based statistics  */
    extern int output_file_num;
#endif


/**
 * @class TallyBank TallyBank.h "pinspec/src/TallyBank.h"
 * @brief The TallyBank contains all tallies for a simulation.
 * @details The TallyBank ensures accurate error checking for tallies
 *          used in a simulation. The TallyBank stores tallies in the
 *          appropriate hash table with keys corresponding to the isotope,
 *          material, region and/or geometry in which a tally is to be applied.
 *          The TallyBank can iterate through all of the approrpirate tallies
 *          for a neutron and make a tally in each one, in addition to providing
 *          other functionality to the main Monte Carlo kernel.
 */
class TallyBank {
private:
    /**
     * @brief TallyBank constructor.
     */
    TallyBank() { }

    /**
     * @brief An overloaded assignment function to allow for static
     *        referencing of the TallyBank class.
     * @return a pointer to the static TallyBank class 
     */
    TallyBank &operator=(const TallyBank &) { return *this; }

    /** Container of all registered tallies */
    std::set<Tally*> _all_tallies;
    /** Hash table of all tallies registered for the geometry */
    std::map< Geometry*, std::set<Tally*>* > _geometry_tallies;
    /** Hash table of all tallies registered for a region */
    std::map< Region*, std::set<Tally*>* > _region_tallies;
    /** Hash map of all tallies registered for a material */
    std::map< Material*, std::set<Tally*>* > _material_tallies;
    /** Hash map of all tallies registered for an isotope */
    std::map< Isotope*, std::set<Tally*>* > _isotope_tallies;
 
public:
    /**
     * @brief TallyBank destructor.
     */
    ~TallyBank() { }
 
    /**
     * @brief Gets an instance pointer to this static class.
     * @return Returns a pointer to the static TallyBank class object.
     */ 
    static TallyBank *Get() {
        static TallyBank instance;
	return &instance;
    }

    void registerTally(Tally* tally);
    void registerTally(Tally* tally, Geometry* geometry);	 
    void registerTally(Tally* tally, Region* region);	 
    void registerTally(Tally* tally, Material* material);
    void registerTally(Tally* tally, Isotope* isotope);
    void deregisterTally(Tally* tally);

    void initializeBatchTallies(int num_batches);
    bool isPrecisionTriggered();
    void incrementNumBatches(int num_batches);
    void computeBatchStatistics();
    void computeScaledBatchStatistics(float scale_factor);
    void outputBatchStatistics();
    void tally(neutron* neutron);

    void clearTallies();
};

#endif /* TALLYBANK_H_ */
