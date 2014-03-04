/**
 * @file TallyFactory.h
 * @brief The TallyFactory class.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 19, 2013
 */

#ifndef TALLYFACTORY_H_
#define TALLYFACTORY_H_

#include "Tally.h"


/**
 * @class TallyFactory TallyFactory.h "pinspec/src/TallyFactory.h"
 * @brief A utility class for creating instances of tallies.
 */
class TallyFactory {

private:
    /**
     * @brief TallyFactory constructor.
     */
    TallyFactory() { }

    /**
     * @brief Assignment operator for static referencing of the TallyFactory.
     * @param & the TallyFactory static class object
     * @return a pointer to the TallyFactory static class object
     */
    TallyFactory &operator=(const TallyFactory &) { return *this; }
 
public:
    ~TallyFactory() { }

    /**
     * @brief Returns a static instance of the TallyFactory class.
     * @return a pointer to the static TallyFactory class
     */
    static TallyFactory *Get() {
        static TallyFactory instance;
	return &instance;
    }

    Tally* createTally(Geometry* geometry, tallyType tally_type, 
		       char* tally_name=(char*)""); 
    Tally* createTally(Region* region, tallyType tally_type, 
		       char* tally_name=(char*)"");
    Tally* createTally(Material* material, tallyType tally_type, 
		       char* tally_name=(char*)"");
    Tally* createTally(Isotope* isotope, tallyType tally_type, 
		       char* tally_name=(char*)"");
};

#endif /* TALLYFACTORY_H_ */
