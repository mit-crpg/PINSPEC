from pinspec import *
from numpy import *
from plotter import *
import matplotlib.pyplot as plt


def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 10000
    num_threads = 4
    log_setlevel(INFO)


    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')

    # Plot the microscopic cross sections for each isotope
    plotMicroXS(u235, ['capture', 'elastic', 'fission', 'absorption'])
    plotMicroXS(u238, ['capture', 'elastic', 'fission', 'absorption'])
    plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    
    
    # Define materials
    mix = Material()
    mix.setMaterialName('fuel moderator mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(u238, 0.50)
    mix.addIsotope(u235, .025)
    
    log_printf(INFO, "Added isotopes")

    # Plot the mixture macroscopic cross sections
    plotMacroXS(mix, ['capture', 'elastic', 'fission', 'absorption', 'total'])
    
    # Define regions
    region_mix = Region()
    region_mix.setRegionName('infinite medium fuel/moderator mix')
    region_mix.setRegionType(INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, "Made mixture region")
    

	# Define tallies - give them to Regions, Materials, or Isotopes
	# This part is really where we need to know how to pass float
    # arrays to/from SWIG

    # Create a tally for the flux
    flux = Tally('total flux', REGION, FLUX)
    flux.generateBinEdges(1E-7, 1E7, 1000, LOGARITHMIC)
    region_mix.addTally(flux)

    ############################################################################
    #EXAMPLE: How to set tally bin edges 
    ############################################################################
    # Create a tally for the absorption rate
    abs_rate = Tally('absorption rate', MATERIAL, ABSORPTION_RATE)
    abs_rate_bin_edges = array([0.1, 1., 5., 10., 100., 1000.])
    abs_rate.setBinEdges(abs_rate_bin_edges)
    mix.addTally(abs_rate)

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    log_printf(INFO, "Made geometry")

	# Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

    log_printf(INFO, "Ran Monte Carlo")
    
    # plot the flux
    plotFlux(flux)

    # Dump batch statistics to output files to some new directory - gives segmentation fault right now
    # geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')


if __name__ == '__main__':
    
    main()
