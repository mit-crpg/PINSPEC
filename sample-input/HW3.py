# This input file corresponds to 22.211 HW3, where we have: 
# slowing down MC, U238 absorption resonance xs (and 0.1 above the 
# energy range). Scattering xs is asymptotic down-scattering. 
import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter


def main():

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 10000
    num_threads = 4
    log_setlevel(INFO)

    # Call SLBW to create XS
    filename = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
    T=300 #Temp in Kelvin of target nucleus
    SLBW.SLBWXS(filename,T)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
 
    # Define materials
    mix = Material()
    mix.setMaterialName('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc') # we do not have a density from last year.
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(u238, 0.50)
    # mix.addIsotope(u235, .025) # we do not have U235 from last year. 

    log_printf(INFO, 'Added isotopes')

    # Define regions
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, 'Made mixture region')
        
    # plot the fission spectrum the CDF
    plotter.plotFissionSpectrum()

    #Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScatteringPDF(h1)

    # Create a tally for the flux
    flux = Tally('total flux', GEOMETRY, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC)

    ############################################################################
    #EXAMPLE: How to set tally bin edges 
    ############################################################################
    # Create a tally for the absorption rate
    abs_rate = Tally('absorption rate', MATERIAL, MICRO_ABSORPTION_RATE)
    abs_rate_bin_edges = numpy.array([0.1, 1., 5., 10., 100., 1000.])
    abs_rate.setBinEdges(abs_rate_bin_edges)
    mix.addTally(abs_rate)

    # Set a precision trigger: tells simulation to run until maximum relative
    # error is less than the trigger value (4E-3)
    abs_rate.setPrecisionTrigger(RELATIVE_ERROR, 3E-3)

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.addTally(flux)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    log_printf(INFO, 'Made geometry')

    # Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

    
    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')

    # Plotting  
    plotter.plotFlux(flux)
    plotter.plotMicroXS(u235, ['capture', 'absorption'])
    plotter.plotMicroXS(u238, ['capture', 'fission', 'elastic', \
                                   'total', 'absorption'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    plotter.plotMacroXS(mix, ['capture', 'elastic', 'fission', \
                                  'absorption', 'total'])
    

if __name__ == '__main__':

    main()  

