import matplotlib.pyplot as matplt
import numpy
from pinspec import *
import plotter


def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 1000
    num_threads = 8
    log_setlevel(INFO)


    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')

    # Plot the microscopic cross sections for each isotope
    plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission', 'absorption'])
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission', 'absorption'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    
    
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
    plotter.plotMacroXS(mix, ['capture', 'elastic', 'fission', \
                                            'absorption', 'total'])
    
    # Define regions
    region_mix = Region()
    region_mix.setRegionName('infinite medium fuel/moderator mix')
    region_mix.setRegionType(INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, "Made mixture region")

    ############################################################################
    #EXAMPLE: How to plot a fission spectrum CDF 
    ############################################################################
    fissioner = Fissioner()
    fissioner.setNumBins(10000)
    fissioner.setEMax(20)
    fissioner.buildCDF()
    cdf = fissioner.retrieveCDF(fissioner.getNumBins())
    cdf_energies = fissioner.retrieveCDFEnergies(fissioner.getNumBins())
        
    fig = matplt.figure()
    matplt.plot(cdf_energies, cdf)
    matplt.xscale('log')
    matplt.xlabel('Energy [ev]')
    matplt.title('Watt Spectrum CDF')
    matplt.savefig('fission_spectrum_cdf.png')

    ############################################################################
    #EXAMPLE: How to plot a fission spectrum - this is just rough and dirty
    ############################################################################
    num_samples = 100000
    emitted_energies = numpy.zeros(num_samples)
    for i in range(num_samples):
        emitted_energies[i] = fissioner.emitNeutroneV()

    fig = matplt.figure()
    matplt.hist(emitted_energies, 100)
    matplt.savefig('fission_spectrum.png')    

    

	# Define tallies - give them to Regions, Materials, or Isotopes
	# This part is really where we need to know how to pass float
    # arrays to/from SWIG

    # Create a tally for the flux
    flux = Tally('total flux', REGION, FLUX)
    flux.generateBinEdges(1E-3, 2E7, 1000, LOGARITHMIC)
    region_mix.addTally(flux)

    ############################################################################
    #EXAMPLE: How to set tally bin edges 
    ############################################################################
    # Create a tally for the absorption rate
    abs_rate = Tally('absorption rate', MATERIAL, ABSORPTION_RATE)
    abs_rate_bin_edges = numpy.array([0.1, 1., 5., 10., 100., 1000.])
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
    plotter.plotFlux(flux)
    

    # Dump batch statistics to output files to some new directory - gives segmentation fault right now
    # geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')


if __name__ == '__main__':
    
    main()
