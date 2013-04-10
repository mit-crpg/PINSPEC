/*
 * Material.h
 *
 *  Created on: Mar 26, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "Isotope.h"

typedef enum densityUnits {
	GRAM_CM3,
	NUM_CM3
} densityUnit;


class Material {

private:
	char* _material_name;
	float _material_density;
	float _material_number_density;
	float _material_atomic_mass;
	float _buckling_squared;
    float _volume;
    densityUnit _density_unit;

	/* Map of number density and isotope pointers */
	std::map<char*, std::pair<float, Isotope*> > _isotopes;
	std::map<Isotope*, float> _isotopes_AO;

public:
	Material(char* material_name);
	virtual ~Material();
	
	/* getters */
	char* getMaterialName();
	float getMaterialNumberDensity();
	float getDensity();
	float getIsotopeNumDensity(Isotope* isotope);
	bool containsIsotope(Isotope* isotope);
	float getBucklingSquared();
    float getVolume();

    int getNumXSEnergies(char* xs_type) const;

	float getTotalMacroXS(float energy);
	float getTotalMacroXS(int energy_index);
	
	float getElasticMacroXS(float energy);
	float getElasticMacroXS(int energy_index);
	
	float getAbsorptionMacroXS(float energy);
	float getAbsorptionMacroXS(int energy_index);
	
	float getCaptureMacroXS(float energy);
	float getCaptureMacroXS(int energy_index);
	
	float getFissionMacroXS(float energy);
	float getFissionMacroXS(int energy_index);
	
	float getTransportMacroXS(float energy);
	float getTransportMacroXS(int energy_index);

    /* IMPORTANT: The following two class method prototypes must not be changed
     * without changing Geometry.i to allow for the data arrays to be transformed
     * into numpy arrays */
    void retrieveXSEnergies(float* energies, int num_xs, 
                                    char* xs_type) const;
    void retrieveXS(float* xs, int num_xs, char* xs_type);
		
	/* setters */
	void setDensity(float density, const char* unit);
	void setAtomicMass(float atomic_mass);
	void setBucklingSquared(float buckling_squared);
    void incrementVolume(float volume);
	void addIsotope(Isotope *isotope, float atomic_ratio);
	Material *clone();

	void sampleIsotope(neutron* neutron);
	void collideNeutron(neutron* neutron);
};


#endif /* MATERIAL_H_ */
