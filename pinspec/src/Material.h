/**
 * @file Material.h
 * @brief The Material class.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 26, 2012
 *
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "Isotope.h"

/**
 * @enum densityUnits
 * @brief Units in which density for isotopes may be expressed
 */

/**
 * @var densityUnit
 * @brief A unit in which density for isotopes may be expressed
 */
typedef enum densityUnits {
	GRAM_CM3,
	NUM_CM3
} densityUnit;



/**
 * @class Material Material.h "pinspec/src/Material.h"
 * @brief The Material class represents a collection of isotope objects.
 * @details The Material class represents a collection of isotope objects
 *          and samples isotopes for collisions with neutrons.
 */
class Material {

private:
    /** The name of the material arbitrarily defined by the user*/
    char* _material_name;
    /** The material density in g/cc */
    float _material_density;
    /** The material number density in num/cc */
    float _material_number_density;
    /** The material atomic mass (sum of all isotope masses according to
     * their abundance in the material */
    float _material_atomic_mass;

    /** The geometric buckling squared - used to sample leakage */
    float _buckling_squared;
    /** The total volume of all regions containing this material */
    float _volume;

    /** Map relating isotope name to number density / isotope pairs */
    std::map<char*, std::pair<float, Isotope*> > _isotopes;
    /** Map relating isotope pointers to atomic ratios within the material */
    std::map<Isotope*, float> _isotopes_AO;

    /** The units for the material's density (ie, 'g/cc' or 'at/cc') */
    densityUnit _density_unit;

public:
    Material(char* material_name);
    virtual ~Material();
	
    /* getters */
    char* getMaterialName();
    float getMaterialNumberDensity();
    Isotope* getIsotope(char* isotope);
    float getDensity();
    float getIsotopeNumDensity(char* isotope);
    bool containsIsotope(Isotope* isotope);
    float getBucklingSquared();
    float getVolume();

    int getNumXSEnergies(char* xs_type);

    float getTotalMacroXS(float energy);
    float getTotalMacroXS(int energy_index);
    float getTotalMicroXS(float energy);
    float getTotalMicroXS(int energy_index);
	
    float getElasticMacroXS(float energy);
    float getElasticMacroXS(int energy_index);
    float getElasticMicroXS(float energy);
    float getElasticMicroXS(int energy_index);

    float getAbsorptionMacroXS(float energy);
    float getAbsorptionMacroXS(int energy_index);
    float getAbsorptionMicroXS(float energy);
    float getAbsorptionMicroXS(int energy_index);
	
    float getCaptureMacroXS(float energy);
    float getCaptureMacroXS(int energy_index);
    float getCaptureMicroXS(float energy);
    float getCaptureMicroXS(int energy_index);
	
    float getFissionMacroXS(float energy);
    float getFissionMacroXS(int energy_index);
    float getFissionMicroXS(float energy);
    float getFissionMicroXS(int energy_index);
	
    float getTransportMicroXS(float energy);
    float getTransportMicroXS(int energy_index);
    float getTransportMacroXS(float energy);
    float getTransportMacroXS(int energy_index);

    /* IMPORTANT: The following two class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be 
     * transformed into numpy arrays */
    void retrieveXSEnergies(float* energies, int num_xs, char* xs_type);
    void retrieveXS(float* xs, int num_xs, char* xs_type);
		
    /* setters */
    void setMaterialName(char* name);
    void setDensity(float density, char* unit);
    void setNumberDensity(float number_density, const char* unit);
    void setAtomicMass(float atomic_mass);
    void setBucklingSquared(float buckling_squared);
    void incrementVolume(float volume);
    void addIsotope(Isotope *isotope, float atomic_ratio);
    
    Material *clone();

    void sampleIsotope(neutron* neutron);
    void collideNeutron(neutron* neutron);
};


#endif /* MATERIAL_H_ */
