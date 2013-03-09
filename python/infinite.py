import matplotlib.pyplot as plt
from pinspec import *
from numpy import *
from scipy import *
#from plotter import *


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

    # EXAMPLE: How to retrieve micro xs from C++ and plot them
    # First, retrieve the xs and x
    u235_num_xs = u235.getNumXSEnergies()
    u235_capture_xs = u235.retrieveXS(u235_num_xs, 'capture')
    u235_elastic_xs = u235.retrieveXS(u235_num_xs, 'elastic')
    u235_fission_xs = u235.retrieveXS(u235_num_xs, 'fission')
    u235_absorption_xs = u235.retrieveXS(u235_num_xs, 'absorption')
    u235_total_xs = u235.retrieveXS(u235_num_xs, 'total')
    #Note: after rescaling, all xs types have the same energy grid so only
    #need to retrieve it for one xs type
    u235_xs_energies = u235.retrieveXSEnergies(u235_num_xs, 'capture')    

    # Plot all of the xs on the same scale
    fig = plt.figure()
    plt.plot(u235_xs_energies, u235_capture_xs, lw=1)
    plt.plot(u235_xs_energies, u235_elastic_xs, lw=1)
    plt.plot(u235_xs_energies, u235_fission_xs, lw=1)
    plt.plot(u235_xs_energies, u235_absorption_xs, lw=1)
    plt.plot(u235_xs_energies, u235_total_xs, lw=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [ev]')
    plt.ylabel('Micro XS [barns]')
    plt.title(h1.getIsotopeType() + ' Micro XS')
    plt.legend(['capture', 'elastic', 'fission', 'absorption', 'total'])
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

    log_printf(INFO, "Added isotopes")
    
    # Define regions
    region_mix = Region()
    region_mix.setRegionName('infinite medium fuel/moderator mix')
    region_mix.setRegionType(INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, "Made mixture region")
    
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

    log_printf(INFO, "Made geometry")

	# Run Monte Carlo simulation
    # geometry.runMonteCarloSimulation();
        
	# Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')
        
	# Plot data

	# Cleanup data - I think you need to do this to avoid a 
	# segmentation fault, but I'm not sure

if __name__ == '__main__':
    
    main()
