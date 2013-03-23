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

	/* Map of number density and isotope pointers */
	std::map<char*, std::pair<float, Isotope*> > _isotopes;
	std::map<Isotope*, float> _isotopes_AO;

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

    int getNumXSEnergies(char* xs_type) const;

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
     * without changing Geometry.i to allow for the data arrays to be transformed
     * into numpy arrays */
    void retrieveXSEnergies(float* energies, int num_xs, 
                                    char* xs_type) const;
    void retrieveXS(float* xs, int num_xs, char* xs_type);
		
	/* setters */
	void setMaterialName(char* name);
	void setDensity(float density, char* unit);
	void setNumberDensity(float number_density);
	void setAtomicMass(float atomic_mass);
	void setBucklingSquared(float buckling_squared);
	void addIsotope(Isotope *isotope, float atomic_ratio);
	Material *clone();

	void sampleIsotope(neutron* neutron);
	void collideNeutron(neutron* neutron);
};


#endif /* MATERIAL_H_ */
