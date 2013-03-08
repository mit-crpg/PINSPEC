import matplotlib.pyplot as plt
from pinspec import *

def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material
    
    log_setlevel(DEBUG)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    
    # Define materials
    mix = Material()
    mix.setMaterialName('fuel moderator mix')
    mix.setDensity(1.0, 'g/cc')
    mix.addIsotope(h1, 2.0)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(u235, 0.03)
    mix.addIsotope(u238, 0.97)

    print 'added isotopes'
   
    # Define regions
    region_mix = Region()
    region_mix.setRegionName('infinite medium fuel/moderator mix')
    region_mix.setRegionType(INFINITE)
    region_mix.setMaterial(mix)

    print 'made region'

	# Define tallies - give them to Regions, Materials, or Isotopes
	# This part is really where we need to know how to pass float
    # arrays to/from SWIG

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.setNeutronsPerBatch(10000)
    geometry.setNumBatches(10)
    geometry.setNumThreads(1)

    print 'made geometry'
        
	# Run Monte Carlo simulation
#    geometry.runMonteCarloSimulation();
        
	# Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')
        
	# Plot data

	# Cleanup data - I think you need to do this to avoid a 
	# segmentation fault, but I'm not sure

if __name__ == '__main__':
    
    main()
