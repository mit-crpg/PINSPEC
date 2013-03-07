import matplotlib.pyplot as plt
from pinspec import *
from numpy import *
from scipy import *


def main():


    # NOTE: If a user is going to homogenize cladding/gap with the moderator
	#       then they must figure out the atom ratio using geometric parameters
	#       and use that when loading the isotope in a material

    log_setlevel('DEBUG')
    
    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')

    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    
    # Define materials
    moderator = Material()
    print 'made material'
    moderator.setMaterialName('moderator')
    print 'set moderator name'
    moderator.setDensity(1.0, 'g/cc')
    print 'set moderator density'
    moderator.addIsotope(h1, 2.0)
    print 'added h1 to moderator'
    moderator.addIsotope(o16, 1.0)

    print 'added isotopes to moderator'
    
    fuel = Material()
    print 'made fuel'
    fuel.setMaterialName('fuel')
    print 'set fuel name'
    fuel.setDensity(10.0, 'g/cc')
    print 'set fuel density'
    fuel.addIsotope(u235, 0.03)
    print 'added u235 to fuel'
    fuel.addIsotope(u238, 0.97)
    print 'added u238 to fuel'
    fuel.addIsotope(o16, 2.0)
    print 'added 016 to fuel'
    
    print 'added isotopes to fuel'
   
    # Define regions
    region_mod = Region()
    region_mod.setRegionType(MODERATOR)
    region_mod.setMaterial(moderator)
    region_mod.setRegionName('moderator')
    
    print 'set moderator'

    region_fuel = Region()
    region_fuel.setRegionType(FUEL)
    region_fuel.setMaterial(fuel)
    region_fuel.setRegionName('fuel')

    print 'set fuel'
    
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
    geometry.runMonteCarloSimulation();

	# Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('DirectoryName')

	# Plot data

	# Cleanup data - I think you need to do this to avoid a 
	# segmentation fault, but I'm not sure
    del geometry
	

if __name__ == '__main__':
    
    main()

