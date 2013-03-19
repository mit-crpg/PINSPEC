import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

def main():

    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 1000
    num_threads = 4
    
    log_setlevel(INFO)

    setXSLibDirectory('../xs-lib/')   # This is also a default, but set it as example

    py_printf('INFO', 'Creating SLBW xs')    
    
    # Call SLBW to create XS
    filename = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
    T=300 #Temp in Kelvin of target nucleus
    SLBW.SLBWXS(filename,T)

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

    zr90 = Isotope('Zr-90')
 
    # Define materials
    mix = Material()
    mix.setMaterialName('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(u238, 0.50)
    mix.addIsotope(u235, .025)
    mix.addIsotope(zr90, .16)
    
    py_printf('INFO', 'Added isotopes')

    plotter.plotMacroXS(mix, ['total', 'capture', 'elastic', 'absorption'])
    
    # Define regions
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    py_printf('INFO', 'Made mixture region')
        
    # plot the fission spectrum the CDF
    #plotter.plotFissionSpectrum()

    #Plot the thermal scattering kernel PDFs and CDFs
    #plotter.plotThermalScattering(h1)

    # Create a tally for the flux
    flux = Tally('total flux', GEOMETRY, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC)

    ############################################################################
    #EXAMPLE: How to set tally bin edges 
    ############################################################################
    # Create a tally for the RI
    abs_rate = Tally('micro absorption rate', MATERIAL, MICRO_ABSORPTION_RATE)
    flux_RI = Tally('flux RI', MATERIAL, FLUX)
    abs_rate_bin_edges = numpy.array([0.1, 1., 6., 10., 25., 50., 100.])
    abs_rate.setBinEdges(abs_rate_bin_edges)
    flux_RI.setBinEdges(abs_rate_bin_edges)
    mix.addTally(abs_rate)
    mix.addTally(flux_RI)

    # Set a precision trigger: tells simulation to run until maximum relative
    # error is less than the trigger value (4E-3)
    abs_rate.setPrecisionTrigger(RELATIVE_ERROR, 7E-3)

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.addTally(flux)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    py_printf('INFO', 'Made geometry')

	# Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')

    # Plot fluxes
    plotter.plotFluxes([flux])
    plotter.plotFluxes([flux])
    plotter.plotFluxes([flux])
    
    # Compute the resonance integrals
    RI = process.RI(flux_RI, abs_rate)
    RI.printRI()
    plotter.plotRI(RI)
    
    # Compute the group cross sections
    groupXS = process.groupXS(flux_RI, abs_rate)
    groupXS.printGroupXS()
    plotter.plotGroupXS(groupXS)

    # Dump batch statistics to output files to some new directory - gives segmentation fault right now
    # geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')


if __name__ == '__main__':

    main()  

