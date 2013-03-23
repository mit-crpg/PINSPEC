/*
 * TallyFactory.h
 *
 *  Created on: Mar 19, 2013
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TALLYFACTORY_H_
#define TALLYFACTORY_H_

#ifdef __cplusplus

#include "Tally.h"


/* Factory for creating instances of Tallies */
class TallyFactory {
private:
	TallyFactory() { }
	TallyFactory(const TallyFactory &) { }
	TallyFactory &operator=(const TallyFactory &) { return *this; }
 
public:
	~TallyFactory() { }
 
	static TallyFactory *Get() {
		static TallyFactory instance;
		return &instance;
	}

	Tally* createTally(char* tally_name, Geometry* geometry, tallyType tally_type); 
	Tally* createTally(char* tally_name, Region* region, tallyType tally_type);
	Tally* createTally(char* tally_name, Material* material, tallyType tally_type);
	Tally* createTally(char* tally_name, Isotope* isotope, tallyType tally_type);
};


#endif

#endif /* TALLYFACTORY_H_ */
