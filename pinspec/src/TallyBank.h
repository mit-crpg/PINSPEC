/*
 * TallyBank.h
 *
 *  Created on: Mar 20, 2013
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TALLYBANK_H_
#define TALLYBANK_H_

#ifdef __cplusplus
#include <set>
#include <map>
#include <utility>
#include "Tally.h"
#include "Region.h"
#include "Material.h"
#include "Isotope.h"
#include "Geometry.h"


/* Factory for creating instances of Tallies */
class TallyBank {
private:
	TallyBank() { }
	TallyBank(const TallyBank &) { }
	TallyBank &operator=(const TallyBank &) { return *this; }

	std::set<Tally*> _all_tallies;
	std::map< Geometry*, std::set<Tally*>* > _geometry_tallies;
	std::map< Region*, std::set<Tally*>* > _region_tallies;
	std::map< Material*, std::set<Tally*>* > _material_tallies;
	std::map< Isotope*, std::set<Tally*>* > _isotope_tallies;
 
public:
	~TallyBank() { }
 
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

#endif

#endif /* TALLYBANK_H_ */
