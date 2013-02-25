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

#include "log.h"
#include "arraycreator.h"
#include "Isotope.h"

typedef enum energyGridType {
	EQUAL,
	LOGARITHMIC,
	OTHER
} energyGridType;


class Material {

private:
	char* _material_name;

	/* Map of number density and isotope pointers */
	std::map<char*, std::pair<float, Isotope*> > _isotopes;
	float _tot_num_density;

	/* Values related to rescaled cross-sections on a uniform energy grid */
	bool _rescaled;
	energyGridType _scale_type;
	int _num_energies;
	float _start_energy;
	float _end_energy;
	float _delta_energy;

public:
	Material();
	virtual ~Material();
	char* getMaterialName();
    float getTotalNumberDensity();
    Isotope* getIsotope(char* isotope);
    float getIsotopeNumDensity(char* isotope);

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

    bool isRescaled();
    int getEnergyGridIndex(float energy);

    void setMaterialName(char* name);
    void addIsotope(Isotope *material, float num_density);

    void rescaleCrossSections(float start_energy, float end_energy,
    						int num_energies, energyGridType scale_type);
    Isotope* sampleIsotope(float energy);
};

#endif /* MATERIAL_H_ */
