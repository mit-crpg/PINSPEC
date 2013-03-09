import numpy as np
import matplotlib.pyplot as plt
from pinspec import *
import time

def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material
    
    log_setlevel(INFO)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')

#    time.sleep(3)

    # EXAMPLE: How to retrieve micro xs from C++ and plot them
    # First, retrieve the xs and x
    h1_num_xs = h1.getNumXSEnergies()
    h1_capture_xs = h1.retrieveXS(h1_num_xs, 'capture')
    h1_elastic_xs = h1.retrieveXS(h1_num_xs, 'elastic')
    h1_fission_xs = h1.retrieveXS(h1_num_xs, 'fission')
    h1_absorption_xs = h1.retrieveXS(h1_num_xs, 'absorption')
    h1_total_xs = h1.retrieveXS(h1_num_xs, 'total')
    #Note: after rescaling, all xs types have the same energy grid so only
    #need to retrieve it for one xs type
    h1_xs_energies = h1.retrieveXSEnergies(h1_num_xs, 'capture')
    print 'h1_num_xs = ' + str(h1_num_xs)
    print 'h1_capture_xs = ' + str(h1_capture_xs)

    fig = plt.figure()
    plt.plot(h1_xs_energies, h1_elastic_xs, lw=3)
    plt.xscale('log')
#    plt.ylim([h1_xs_energies[0], h1_xs_energies[h1_xs_energies.size-1]])
    plt.xlabel('Energy [ev]')
    plt.ylabel('Micro XS [barns]')
    plt.title(h1.getIsotopeType() + ' Micro XS')
    fig.show()
    fig.savefig('h1_capture_xs.png')

    
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
