/**
 * @file Neutron.h
 * @brief The neutron C structure 
 * @author William Boyd (wboyd@mit.edu)
 * @date March 13, 2012
 *
 */


#ifndef NEUTRON_H_
#define NEUTRON_H_

class Region;	
class Material;
class Isotope;


/**
 * @struct neutron
 * @brief Represents a neutron in a PINSPEC simulation.
 * @details The neutron struct is the fundamental unit of data needed to
 *          represent a neutron in a PINSPEC simulation. This struct is useful
 *          for passing a neutron from object to object, such as from the
 *          Geometry to a Region to a Material and to an Isotope to update
 *          the neutron's energy from a collision. 
 */
struct neutron {
    /** this neutron's batch number in the simulation */
    int _batch_num;
    /** the x-coordinate of this neutron's location */
    float _x;
    /** the y-coordinate of this neutron's location */
    float _y;
    /** the z-coordinate of this neutron's location */
    float _z;
    /** the component of this neutron's velocity unit vector along the x-axis */
    float _u;
    /** the component of this neutron's velocity unit vector along the y-axis */
    float _v;
    /** the component of this neutron's velocity unit vector along the z-axis */
    float _w;

    /** the cosine of the polar angle \f$\theta\f$ of this neutron's direction 
     * vector: \f$\mu = cos(\theta)\f$
     */
    float _mu;

    /** the azimuthal angle \f$\phi\f$ of this neutron's direction vector in 
     * the xy-plane */
    float _phi;

    /** the neutron's energy in eV following its most recent collision */
    float _energy;

    /** a boolean representing whether this neutron has been absorbed or not */
    bool _alive;

    /** a pointer to the Regon in which the neutron most recently collided */
    Region* _region;

    /** a pointer to the Material with which the neutron most recently
      * collided 
      */
    Material* _material;

    /** a pointer to the Isotope with which the neutron most recently 
      * collided */
    Isotope* _isotope;

    /** the neutron's energy in eV prior to its most recent collision */
    float _old_energy;

    /** the total macroscopic cross-section in \f$cm^{-1}\f$ of the material 
     * in which the neutron collided at its collision energy */
    float _total_xs;
};

neutron* initializeNewNeutron();


#endif /* NEUTRON_H_ */
