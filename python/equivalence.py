import matplotlib.pyplot as plt
from pinspec import *
from numpy import *
from scipy import *
from plotter import *


def main():


    # NOTE: If a user is going to homogenize cladding/gap with the moderator
	#       then they must figure out the atom ratio using geometric parameters
	#       and use that when loading the isotope in a material

    log_setlevel(INFO)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    
    # Define moderator material
    moderator = Material()
    moderator.setMaterialName('moderator')
    moderator.setDensity(1.0, 'g/cc')
    moderator.addIsotope(h1, 2.0)
    moderator.addIsotope(o16, 1.0)

    log_printf(INFO, "Added isotopes to moderator")

    # Define fuel material
    fuel = Material()
    fuel.setMaterialName('fuel')
    fuel.setDensity(10.0, 'g/cc')
    fuel.addIsotope(u235, 0.03)
    fuel.addIsotope(u238, 0.97)
    fuel.addIsotope(o16, 2.0)
    
    log_printf(INFO, "Added isotopes to fuel")
    
    # Define moderator region
    region_mod = Region()
    region_mod.setRegionName('moderator')
    region_mod.setRegionType(MODERATOR)
    region_mod.setMaterial(moderator)
    
    log_printf(INFO, "Made moderator region")
    
    # Define fuel region
    region_fuel = Region()
    region_fuel.setRegionName('fuel')
    region_fuel.setRegionType(FUEL)
    region_fuel.setMaterial(fuel)

    log_printf(INFO, "Made fuel region")
    
	# Define tallies - give them to Regions, Materials, or Isotopes
	# This part is really where we need to know how to pass float
    # arrays to/from SWIG

	# Two region homogeneous equivalence parameters
    radius_fuel = 0.4096;
    dancoff = 0.277;
    sigma_e = 1.0 / (2.0*radius_fuel);
    A = (1.0 - dancoff) / dancoff;
    alpha1 = ((5.0*A + 6.0) - sqrt(A*A + 36.0*A + 36.0)) / (2.0*(A+1.0));
    alpha2 = ((5.0*A + 6.0) + sqrt(A*A + 36.0*A + 36.0)) / (2.0*(A+1.0));
    beta = (((4.0*A + 6.0) / (A + 1.0)) - alpha1) / (alpha2 - alpha1);

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(HOMOGENEOUS_EQUIVALENCE)
    geometry.setTwoRegionPinCellParams(sigma_e, beta, alpha1, alpha2)
    geometry.addRegion(region_mod)
    geometry.addRegion(region_fuel)
    geometry.setNeutronsPerBatch(10000)
    geometry.setNumBatches(10)
    geometry.setNumThreads(1)

	# Run Monte Carlo simulation
    # geometry.runMonteCarloSimulation();

	# Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Equivalence_MC_Statistics', 'test')

	# Plot data

	# Cleanup data - I think you need to do this to avoid a 
	# segmentation fault, but I'm not sure
    del geometry
	

if __name__ == '__main__':
    
    main()

