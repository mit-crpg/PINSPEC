import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

def main():
    
    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 100000
    num_threads = 4
    setOutputDirectory('Equivalence')

    log_setlevel(INFO)

    py_printf('TITLE', 'Simulating an infinite medium homogenized pin-cell')

    py_printf('INFO', 'Initializing isotopes...')

    # Define isotopes
    h1 = Isotope('H-1')
    b10 = Isotope('B-10')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    zr90 = Isotope('Zr-90')    

    py_printf('INFO', 'Initializing fuel-moderator mix material...')

    # Define materials
    mix = Material('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(b10, .0000001)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(u238, 0.01)
    mix.addIsotope(u235, .0025)
    mix.addIsotope(zr90, .16)

    py_printf('INFO', 'Initializing fuel-moderator mix region...')
    
    # Define region
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    py_printf('INFO', 'Initializing the geometry...')

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    py_printf('INFO', 'Initializing flux tally...')

    # Create a tally for the flux
    flux = TallyFactory.createTally(region_mix, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC) 
    
    # Set a precision trigger: tells simulation to run until maximum relative
    # error is less than the trigger value (2E-2)
    flux.setPrecisionTrigger(RELATIVE_ERROR, 2E-2)

	# Register the tallies
    TallyBank.registerTally(flux)

    # Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

    py_printf('INFO', 'Plotting microscopic and macroscopic cross-sections...')

    # Plot the microscopic cross sections for each isotope
    plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission'])
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    plotter.plotMacroXS(mix, ['capture', 'elastic', 'fission'])

    py_printf('INFO', 'Plotting fluxes...')

    # Plot flux
    plotter.plotFlux(flux)

    py_printf('INFO', 'Writing tally batch statistics to output file...')

    # Dump batch statistics to output files to some new directory
    TallyBank.outputBatchStatistics()

    py_printf('TITLE', 'Finished')


if __name__ == '__main__':

    main()  

