#include "TallyBank.h"

#define TALLYBANK_C

int output_file_num = 0;


/**
 * @brief Register a Tally for the Tally's domain.
 * @details Registers a tally to be used for any collision within the 
 *          Tally's domain (ie, ISOTOPE, MATERIAL, REGION, GEOMETRY).
 * @param tally the tally to register
 */
void TallyBank::registerTally(Tally* tally) {

    /* Track tallies with a GEOMETRY domain in geometry */
    if (tally->getTallyDomainType() == GEOMETRY) {
        GeometryTally* geometry_tally = static_cast<GeometryTally*>(tally);
        registerTally(geometry_tally, geometry_tally->getGeometry());
    }

    /* Track tallies with a REGION domain in region */
    else if (tally->getTallyDomainType() == REGION) {
        RegionTally* region_tally = static_cast<RegionTally*>(tally);
	registerTally(region_tally, region_tally->getRegion());
    }

    /* Track tallies with a MATERIAL domain in material */
    else if (tally->getTallyDomainType() == MATERIAL) {
        MaterialTally* material_tally = static_cast<MaterialTally*>(tally);
	registerTally(material_tally, material_tally->getMaterial());
    }

    /* Track tallies with a ISOTOPE domain in isotope */
    else if (tally->getTallyDomainType() == ISOTOPE) {
        IsotopeTally* isotope_tally = static_cast<IsotopeTally*>(tally);
	registerTally(isotope_tally, isotope_tally->getIsotope());
    }

    /* Don't track tallies with UNDEFINED domain (only for DERIVED tallies) */
    else if (tally->getTallyDomainType() == UNDEFINED)
        log_printf(ERROR, "Unable to register DERIVED type tally %s with the"
                    " TallyBank", tally->getTallyName());	
}


/**
 * @brief Register a Tally for the geometry.
 * @details Registers a tally to be used whenever a collision occurs
 *          anywhere within the geometry.
 * @param tally the tally to register
 * @param geometry a pointer to geometry object within which we should tally
 */
void TallyBank::registerTally(Tally* tally, Geometry* geometry) {

    /* Track tallies with a GEOMETRY domain in geometry */
    if (tally->getTallyDomainType() == GEOMETRY) {

        /* Add this tally to the geometry's tally registry */
        /* If this geometry is not registered yet, create a new pair for it*/
        if (_geometry_tallies.find(geometry) == _geometry_tallies.end()) {
            std::set<Tally*>* tally_set = new std::set<Tally*>;
	    tally_set->insert(tally);
	    _geometry_tallies[geometry] = tally_set;
        }

        /* Register this tally for an existing geometry in the TallyBank */
        else
            _geometry_tallies.find(geometry)->second->insert(tally);

        /* Add the tally to the global registry */
        _all_tallies.insert(tally);

        log_printf(INFO, "Registered tally %s with the TallyBank for"
                        " the geometry", tally->getTallyName());
        }


        /* Track tallies with a REGION domain in the geometry */
        else if (tally->getTallyDomainType() == REGION) {

        /* Type cast this tally as a region tally */
        RegionTally* region_tally = static_cast<RegionTally*>(tally);
        registerTally(region_tally, region_tally->getRegion());
    }


    /* Track tallies with a MATERIAL domain in the geometry */
    else if (tally->getTallyDomainType() == MATERIAL) {

        /* Type cast this tally as a material tally */
        MaterialTally* material_tally = static_cast<MaterialTally*>(tally);
        registerTally(material_tally, material_tally->getMaterial());
    }


    /* Track tallies with an ISOTOPE domain in the geometry */
    else if (tally->getTallyDomainType() == ISOTOPE) {
	
        /* Type cast this tally as an isotope tally */
	IsotopeTally* isotope_tally = static_cast<IsotopeTally*>(tally);
	registerTally(isotope_tally, isotope_tally->getIsotope());
    }

    /* Don't track tallies with UNDEFINED domain (only for DERIVED tallies) */
    else if (tally->getTallyDomainType() == UNDEFINED)
        log_printf(ERROR, "Unable to register DERIVED type tally %s with the"
                        " TallyBank", tally->getTallyName());

}


/**
 * @brief Register a Tally for a region.
 * @details Registers a tally to be used whenever a collision occurs
 *          anywhere within the user-specified region.
 * @param tally the tally to register
 * @param region a pointer to the region within which we should tally
 */
void TallyBank::registerTally(Tally* tally, Region* region) {

    /* We are unable to track tallies with a GEOMETRY domain for a region */
    if (tally->getTallyDomainType() == GEOMETRY)
	log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" for a region since it is a GEOMETRY "
			"type Tally", tally->getTallyName());

    /* Track tallies with a REGION domain in this region */
    else if (tally->getTallyDomainType() == REGION) {

        /* Type cast this tally as a region tally */
        RegionTally* region_tally = static_cast<RegionTally*>(tally);

        /* Don't register tally in region if tally's region doesn't match */
        if (region_tally->getRegion() != region)
             log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in region %s since it is a REGION type tally "
			" for region %s", region_tally->getTallyName(), 
			region->getName(), 
			region_tally->getRegion()->getName());
    }

    /* Track tallies with a MATERIAL domain in this region */
    else if (tally->getTallyDomainType() == MATERIAL) {

    /* Type cast this tally as a material tally */
    MaterialTally* material_tally = static_cast<MaterialTally*>(tally);

    /* Don't register tally if tally's material doesn't match */
    if (material_tally->getMaterial() != region->getMaterial())
        log_printf(ERROR, "The TallyBank is unable to register tally %s"
		   " in region %s with material %s since it is a "
		   "MATERIAL type tally for material %s", 
		   material_tally->getTallyName(), 
		   region->getName(), 
		   region->getMaterial()->getMaterialName(),
		   material_tally->getMaterial()->getMaterialName());
    }


    /* Track tallies with an ISOTOPE domain in this region */
    else if (tally->getTallyDomainType() == ISOTOPE) {
	
    /* Type cast this tally as a material tally */
    IsotopeTally* isotope_tally = static_cast<IsotopeTally*>(tally);

    /* Don't register tally in region if tally's isotope isn't 
     * contained by the region's material */
    if (!region->containsIsotope(isotope_tally->getIsotope()))
        log_printf(ERROR, "The TallyBank is unable to register tally %s"
		   " in region %s since it is an ISOTOPE type tally for"
		   " isotope %s which is not contained in material %s", 
		   isotope_tally->getTallyName(), 
		   region->getName(), 
		   isotope_tally->getIsotope()->getIsotopeName(),
		   region->getMaterial()->getMaterialName());
    }

    /* Don't track tallies with UNDEFINED domain (only for DERIVED tallies) */
    else if (tally->getTallyDomainType() == UNDEFINED)
        log_printf(ERROR, "Unable to register DERIVED type tally %s with the"
                        " TallyBank", tally->getTallyName());

    /* Don't track tallies of UNDEFINED domain (only for DERIVED tallies) */
    else if (tally->getTallyDomainType() == UNDEFINED)
        log_printf(ERROR, "Unable to register DERIVED type tally %s with"
		       " the TallyBank", tally->getTallyName());

    /* If the tally and region passed all tests, register the region
     * and add this tally to the region's tally registry */

    /* If this region is not registered yet, create a new pair for it*/
    if (_region_tallies.find(region) == _region_tallies.end()) {
        std::set<Tally*>* tally_set = new std::set<Tally*>;
        tally_set->insert(tally);
        _region_tallies[region] = tally_set;
    }

    /* Register this tally for an existing Region in the TallyBank */
    else
        _region_tallies.find(region)->second->insert(tally);

    /* Add the tally to the global registry */
    _all_tallies.insert(tally);

    log_printf(INFO, "Registered tally %s with the TallyBank for region %s", 		        tally->getTallyName(), region->getName());
}


/**
 * @brief Register a Tally for a material.
 * @details Registers a tally to be used whenever a collision occurs
 *          anywhere within the material.
 * @param tally the tally to register
 * @param material a pointer to the material object within which we should tally
 */
void TallyBank::registerTally(Tally* tally, Material* material) {

    /* We are unable to track tallies with a GEOMETRY domain for a material */
    if (tally->getTallyDomainType() == GEOMETRY)
        log_printf(ERROR, "The TallyBank is unable to register tally %s"
		   " for a material since it is a GEOMETRY "
		   "type tally", tally->getTallyName());

    /* Track tallies with a REGION domain in this material */
    else if (tally->getTallyDomainType() == REGION)
        log_printf(ERROR, "The TallyBank is unable to register tally %s"
		   " for a material since it is a REGION "
		   "type tally", tally->getTallyName());

    /* Track tallies with a MATERIAL domain in this material */
    else if (tally->getTallyDomainType() == MATERIAL) {

        /* Type cast this tally as a material tally */
        MaterialTally* material_tally = static_cast<MaterialTally*>(tally);

        /* Don't register tally if tally's material doesn't match */
        if (material_tally->getMaterial() != material)
            log_printf(ERROR, "The TallyBank is unable to register tally %s"
		       " in material %s since it is a MATERIAL type tally "
		       " for material %s", material_tally->getTallyName(), 
		       material->getMaterialName(), 
		       material_tally->getMaterial()->getMaterialName());
    }


    /* Track tallies with an ISOTOPE domain in this region */
    else if (tally->getTallyDomainType() == ISOTOPE) {
	
	/* Type cast this tally as a material tally */
	IsotopeTally* isotope_tally = static_cast<IsotopeTally*>(tally);

	/* Don't register tally in material if tally's isotope isn't 
         * contained by the material */
	if (!material->containsIsotope(isotope_tally->getIsotope()))
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in material %s since it is an ISOTOPE type tally for "
			" isotope %s which is not contained in material %s", 
			isotope_tally->getTallyName(), 
		        material->getMaterialName(), 
			isotope_tally->getIsotope()->getIsotopeName(),
			material->getMaterialName());
	}


        /* If the tally and material passed all tests, register the material
	 * and add this tally to the material's tally registry.
	 * If this material is not registered yet, create a new pair for it*/
        if (_material_tallies.find(material) == _material_tallies.end()) {
	    std::set<Tally*>* tally_set = new std::set<Tally*>;
	    tally_set->insert(tally);
	    _material_tallies[material] = tally_set;
	}

        /* Register this tally for an existing Material in the TallyBank */
	else
	    _material_tallies.find(material)->second->insert(tally);

        /* Add the tally to the global registry */
	_all_tallies.insert(tally);

	log_printf(INFO, "Registered tally %s with the TallyBank for material"
		   " %s", tally->getTallyName(), material->getMaterialName());
}



/**
 * @brief Register a Tally for an isotope.
 * @details Registers a tally to be used whenever a collision occurs
 *          anywhere within the geometry but using only the isotope's
 *          microscopic cross-section data.
 * @param tally the tally to register
 * @param isotope the isotope object within which we should tally
 */
void TallyBank::registerTally(Tally* tally, Isotope* isotope) {

    /* We are unable to track tallies with a GEOMETRY domain for an isotope */
    if (tally->getTallyDomainType() == GEOMETRY)
	log_printf(ERROR, "The TallyBank is unable to register tally %s"
		   " for an isotope since it is a GEOMETRY "
		   "type tally", tally->getTallyName());


    /* We are unable to track tallies with a REGION domain for an isotope */
    else if (tally->getTallyDomainType() == REGION)
        log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" for an isotope since it is a REGION "
			"type tally", tally->getTallyName());

    /* We are unable to track tallies with a MATERIAL domain for an isotope */
    else if (tally->getTallyDomainType() == MATERIAL)
        log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" for an isotope since it is a MATERIAL "
			"type tally", tally->getTallyName());


    /* Track tallies with an ISOTOPE domain in this region */
    else if (tally->getTallyDomainType() == ISOTOPE) {
	
	/* Type cast this tally as a material tally */
	IsotopeTally* isotope_tally = static_cast<IsotopeTally*>(tally);

	/* Don't register tally for isotope if tally's isotope doesnt match */
	if (isotope_tally->getIsotope() != isotope)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in isotope %s since it is an ISOTOPE type tally for "
			" isotope %s", isotope_tally->getTallyName(), 
			isotope->getIsotopeName(), 
			isotope_tally->getIsotope()->getIsotopeName());

	else if (isotope_tally->getTallyType() == FLUX)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in an isotope since it is a FLUX type tally",
			isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == INTERCOLLISION_TIME)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
		       " in an isotope since it is an INTERCOLLISION_TIME"
                        " type tally", isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == LEAKAGE_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in an isotope since it is a LEAKAGE_RATE type tally",
			isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == ELASTIC_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in an isotope since it is a ELASTIC_RATE type tally",
		       isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == GROUP_TO_GROUP_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
		 " in an isotope since it is a GROUP_TO_GROUP_RATE type tally",
		 isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == OUTSCATTER_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
		" in an isotope since it is a OUTSCATTER_RATE type tally",
		       isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == CAPTURE_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
			" in an isotope since it is a CAPTURE_RATE type tally",
			isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == ABSORPTION_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
		" in an isotope since it is a ABSORPTION_RATE type tally",
		isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == FISSION_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
		       " in an isotope since it is a FISSION_RATE type tally",
		       isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == TRANSPORT_RATE)
	    log_printf(ERROR, "The TallyBank is unable to register tally %s"
		" in an isotope since it is a TRANSPORT_RATE type tally", 
		isotope_tally->getTallyName());

	else if (isotope_tally->getTallyType() == DIFFUSION_RATE)
		log_printf(ERROR, "The TallyBank is unable to register tally %s"
			   " in an isotope since it is a DIFFUSION_RATE type "
			   "tally", isotope_tally->getTallyName());
	}

    /* Don't track tallies with UNDEFINED domain (only for DERIVED tallies) */	
    else if (tally->getTallyDomainType() == UNDEFINED)
        log_printf(ERROR, "Unable to register DERIVED type tally %s with the"
                        " TallyBank", tally->getTallyName());

	/* If the tally and isotope passed all tests, register the isotope
	 * and add this tally to the isotope's tally registry */

	/* If this material is not registered yet, create a new pair for it*/
	if (_isotope_tallies.find(isotope) == _isotope_tallies.end()) {
	    std::set<Tally*>* tally_set = new std::set<Tally*>;
	    tally_set->insert(tally);
	    _isotope_tallies[isotope] = tally_set;
	}

	/* Register this tally for an existing isotope in the TallyBank */
	else
	    _isotope_tallies.find(isotope)->second->insert(tally);

	/* Add the tally to the global registry */
	_all_tallies.insert(tally);	
}


/**
 * @brief Un-registers a Tally from the TallyBank.
 * @details Deregisters a tally to be used whenever a collision occurs
 *          anywhere within the tally's domain of geometry, region, material 
 *          or isotope.
 * @param tally the tally to deregister
 */
void TallyBank::deregisterTally(Tally* tally) {

    std::set<Tally*>::iterator set_iter;
    std::map<Geometry*, std::set<Tally*>* >::iterator iter1;
    std::map<Region*, std::set<Tally*>* >::iterator iter2;
    std::map<Material*, std::set<Tally*>* >::iterator iter3;
    std::map<Isotope*, std::set<Tally*>* >::iterator iter4;

    set_iter = _all_tallies.find(tally);
    if (set_iter != _all_tallies.end())
        _all_tallies.erase(set_iter);

    for (iter1 = _geometry_tallies.begin(); iter1 != _geometry_tallies.end(); 
	 ++iter1) {
        set_iter = ((*iter1).second)->find(tally);
        if (set_iter != ((*iter1).second)->end())
            ((*iter1).second)->erase(set_iter);
    }

    for (iter2 = _region_tallies.begin(); iter2 != _region_tallies.end(); 
	 ++iter2) {
        set_iter = ((*iter2).second)->find(tally);
        if (set_iter != ((*iter2).second)->end())
            ((*iter2).second)->erase(set_iter);
    }

    for (iter3 = _material_tallies.begin(); iter3 != _material_tallies.end(); 
	 ++iter3) {
        set_iter = ((*iter3).second)->find(tally);
        if (set_iter != ((*iter3).second)->end())
            ((*iter3).second)->erase(set_iter);
    }

    for (iter4 = _isotope_tallies.begin(); iter4 != _isotope_tallies.end(); 
	 ++iter4) {
        set_iter = ((*iter4).second)->find(tally);
        if (set_iter != ((*iter4).second)->end())
            ((*iter4).second)->erase(set_iter);
    }

    return;
}


/**
 * @brief Checks whether a tally in the TallyBank contains a precision trigger * * @details This method iterates over all tallies which are registered with the 
 *          TallyBank to see if any precision triggers remain.
 * @return Returns true if an active presicion trigger remains
 */
bool TallyBank::isPrecisionTriggered() {

    /* Check if the TallyBank has any tallies with a trigger on precision */
    std::set<Tally*>::iterator iter;

    for (iter = _all_tallies.begin(); iter != _all_tallies.end(); ++iter) {
        if ((*iter)->isPrecisionTriggered())
            return true;
    }

    return false;

}


/**
 * @brief Compute batch statistics for all registered tallies.
 */
void TallyBank::computeBatchStatistics() {

    /* Compute statistics for each of the TallyBank's tallies */
    std::set<Tally*>::iterator iter;

    for (iter = _all_tallies.begin(); iter != _all_tallies.end(); ++iter)
        (*iter)->computeBatchStatistics();

    return;
}


/**
 * @brief Compute scaled batch statistics for all registered tallies.
 * @param scale_factor the value to scaled all of the batch statistics by
 */
void TallyBank::computeScaledBatchStatistics(float scale_factor) {

    std::set<Tally*>::iterator set_iter;
    std::map<Geometry*, std::set<Tally*>* >::iterator iter1;
    std::map<Region*, std::set<Tally*>* >::iterator iter2;
    std::map<Material*, std::set<Tally*>* >::iterator iter3;
    std::map<Isotope*, std::set<Tally*>* >::iterator iter4;

    float volume;

    /* Geometry tallies */
    for (iter1 = _geometry_tallies.begin(); iter1 != _geometry_tallies.end(); 
	 ++iter1) {
        std::set<Tally*> tally_set = (*(*iter1).second);
        volume = (*(*iter1).first).getVolume();

        for (set_iter = tally_set.begin(); set_iter != tally_set.end(); 
	     ++set_iter) {

            if ((*set_iter)->getTallyType() == INTERCOLLISION_TIME)
                (*set_iter)->computeScaledBatchStatistics(scale_factor);
            else
                (*set_iter)->computeScaledBatchStatistics(scale_factor*volume);
        }
    }

    /* Region tallies */
    for (iter2 = _region_tallies.begin(); iter2 != _region_tallies.end(); 
	  ++iter2) {

        std::set<Tally*> tally_set = (*(*iter2).second);
        volume = (*(*iter2).first).getVolume();

        for (set_iter = tally_set.begin(); set_iter != tally_set.end(); 
	     ++set_iter) {

            if ((*set_iter)->getTallyType() == INTERCOLLISION_TIME)
                (*set_iter)->computeScaledBatchStatistics(scale_factor);
            else
                (*set_iter)->computeScaledBatchStatistics(scale_factor*volume);
        }
    }

    /* Material tallies */
    for (iter3 = _material_tallies.begin(); iter3 != _material_tallies.end(); 
	 ++iter3) {

        std::set<Tally*> tally_set = (*(*iter3).second);
        volume = (*(*iter3).first).getVolume();

        for (set_iter = tally_set.begin(); set_iter != tally_set.end(); 
	     ++set_iter) {

            if ((*set_iter)->getTallyType() == INTERCOLLISION_TIME)
                (*set_iter)->computeScaledBatchStatistics(scale_factor);
            else
                (*set_iter)->computeScaledBatchStatistics(scale_factor*volume);
        }
    }

    /* Isotope tallies */
    for (iter4 = _isotope_tallies.begin(); iter4 != _isotope_tallies.end(); 
	 ++iter4) {

        std::set<Tally*> tally_set = (*(*iter4).second);

        for (set_iter = tally_set.begin(); set_iter != tally_set.end(); 
	     ++set_iter)
            (*set_iter)->computeScaledBatchStatistics(scale_factor);
    }

    return;
}


/**
 * @brief Calls each of the Tally class objects in the simulation to output
 *        their tallies and statistics to output files. 
 * @details If a user asks to output the files to a directory which does 
 *          not exist, this method will create the directory.
 */
void TallyBank::outputBatchStatistics() {

    const char* directory = get_output_directory();
    const char* sub_directory = "/tally-statistics";

    std::string full_directory = std::string(directory) + sub_directory;

    /* Check to see if directory exists - if not, create it */
    struct stat st;
    if (!stat(directory, &st) == 0) {
        mkdir(directory, S_IRWXU);
    }

    /* Create tally-statistics subdirectory */
    if (!stat(full_directory.c_str(), &st) == 0) {
        mkdir(full_directory.c_str(), S_IRWXU);
    }

    /* Compute statistics for each of the TallyBank's tallies */
    std::set<Tally*>::iterator iter;

        for (iter = _all_tallies.begin(); iter != _all_tallies.end(); ++iter) {

        std::stringstream filename;

        /* Create filename using the tally's name if it has one */
        if (strlen((*iter)->getTallyName()) != 0)
            filename << (*iter)->getTallyName();

        /* If the tally does not have a name, generate a unique id for it */
        else {
            filename << "tally-" << output_file_num;
            output_file_num++;
        }

        /* Output batch statistics for this tally */
        filename << ".data";
        std::string filename_str = filename.str();
        filename.str("");

        /* Replace spaces with dashes in filename */
        replace(filename_str.begin(), filename_str.end(), ' ', '-');

        /* Make string all lower case */
        std::transform(filename_str.begin(), filename_str.end(), 
                                    filename_str.begin(), ::tolower);

        filename << full_directory << "/" << filename_str;
        (*iter)->outputBatchStatistics(filename.str().c_str());
    }

    return;   
}


/**
 * @brief Tallies a neutron in all appropriate Tally objects.
 * @param neutron the neutron we wish to tally
 */
void TallyBank::tally(neutron* neutron) {

    /* If the neutron has not collided, it was simply transferred between
     * heterogeneous-homogeneous equivalence regions or across a surface
     * in a heterogeneous geometry */
    //FIXME: Tracklength tallies should still account for neutrons which
    //       have not collided
  //    if (neutron->_collided == false)
  //      return;

    std::set<Tally*> tallies;
    std::set<Tally*>::iterator tally_iter;
    std::map<Geometry*, std::set<Tally*>* >::iterator geometry_iter;
    std::map<Region*, std::set<Tally*>* >::iterator region_iter;
    std::map<Material*, std::set<Tally*>* >::iterator material_iter;
    std::map<Isotope*, std::set<Tally*>* >::iterator isotope_iter;

    region_iter = _region_tallies.find(neutron->_region);
    material_iter = _material_tallies.find(neutron->_material);
    isotope_iter = _isotope_tallies.find(neutron->_isotope);

    /* Tally within all tallies registered for the entire geometry */
    for (geometry_iter = _geometry_tallies.begin(); geometry_iter != 
    				_geometry_tallies.end(); ++geometry_iter) {

        tallies = (*(*geometry_iter).second);

        for (tally_iter = tallies.begin(); tally_iter != tallies.end(); 
	     ++tally_iter)
	    (*tally_iter)->tally(neutron);
    }


    /* Tally within all tallies for this neutron's region
     * if tallies are registered for this region */
    if (region_iter != _region_tallies.end()) {
        tallies = *(*region_iter).second;
	for (tally_iter = tallies.begin(); tally_iter !=
					tallies.end(); ++tally_iter) {
   	    (*tally_iter)->tally(neutron);
	}
    }


    /* Tally within all tallies for this neutron's material
    * if tallies are registered for this material */
    if (material_iter != _material_tallies.end()) {
	tallies = *(*material_iter).second;
	for (tally_iter = tallies.begin(); tally_iter !=
						tallies.end(); ++tally_iter)
	    (*tally_iter)->tally(neutron);
	}


    /* Tally within all tallies for this neutron's isotope
     * if tallies are registered for this neutron's isotope */
    if (isotope_iter != _isotope_tallies.end()) {
	tallies = *(*isotope_iter).second;
	for (tally_iter = tallies.begin(); tally_iter !=
						tallies.end(); ++tally_iter)
	    (*tally_iter)->tally(neutron);
	}

}


/**
 * @brief Initializes each registered tally with some number of batches for
 *        batch-based statistics.
 * @param num_batches the total number of batches
 */
void TallyBank::initializeBatchTallies(int num_batches) {

    std::set<Tally*>::iterator iter;

    /* Set the number of batches for all of this Geometry's Tallies */
    for (iter = _all_tallies.begin(); iter != _all_tallies.end(); iter ++)
        (*iter)->setNumBatches(num_batches);

    log_printf(INFO, "TallyBank has initialized %d tallies for %d batches", 
	       _all_tallies.size(), num_batches);
}


/**
 * @brief Increment the total number of batches for batch-based statistics
 *        for each registered Tally object.
 * @param num_batches the number of batches we wish to increment by
 */
void TallyBank::incrementNumBatches(int num_batches) {

    std::set<Tally*>::iterator iter;

    /* Update the number of batches for all of this Geometry's Tallies */
    for (iter = _all_tallies.begin(); iter != _all_tallies.end(); iter ++)
        (*iter)->incrementNumBatches(num_batches);

}


/**
 * @brief Delete all containers (vectors and maps) of Tally objects.
 */
void TallyBank::clearTallies() {

    /* Remove all tally pointers from the TallyBank */
    _all_tallies.clear();
    _geometry_tallies.clear();
    _region_tallies.clear();
    _material_tallies.clear();
    _isotope_tallies.clear();

}
