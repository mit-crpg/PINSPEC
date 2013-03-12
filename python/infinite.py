from pinspec import *
import numpy
import matplotlib.pyplot as matplt    # only need to import this for examples
import plotter
from SLBW import *

def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 100000
    num_threads = 8
    log_setlevel(INFO)

    # Call SLBW to create XS
#    filename = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
#    T=300 #Temp in Kelvin of target nucleus
#    SLBWXS(filename,T)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    
    # Define materials
    mix = Material()
    mix.setMaterialName('fuel moderator mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(u238, 0.50)
    mix.addIsotope(u235, .025)
    
    log_printf(INFO, "Added isotopes")


    # Define regions
    region_mix = Region('infinite medium fuel-moderator mix', INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, "Made mixture region")
        
    # plot the fission spectrum the CDF
    plotter.plotFissionSpectrum()

    #Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScatteringPDF(h1)

    # Create a tally for the flux
    flux = Tally('total flux', REGION, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 2000, LOGARITHMIC)
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

    
    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')

    # Plotting
    plotter.plotFlux(flux)
    plotter.plotMicroXS(u235, ['capture', 'absorption'])
    plotter.plotMicroXS(u238, ['capture', 'absorption'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    plotter.plotMacroXS(mix, ['capture', 'elastic', 'fission', \
                                                    'absorption', 'total'])


if __name__ == '__main__':
    
    main()  
