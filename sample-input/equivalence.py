import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

def main():

    # Set main simulation params
    num_batches = 25
    num_neutrons_per_batch = 1000000
    num_threads = 4
    radius_fuel = 0.4096;
    pitch = 1.26
    dancoff = 0.277;
    setOutputDirectory('Equivalence')

    log_setlevel(INFO)

    py_printf('TITLE', 'Simulating a hetero-homogeneous equivalence pin-cell')

    py_printf('INFO', 'Initializing isotopes...')

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')

    # Neglect thermal scattering in O-16, U-235, U-238
    o16.neglectThermalScattering()
    u235.neglectThermalScattering()
    u238.neglectThermalScattering()

    # Set one group potential elastic scattering xs for u-235
    xs_energies = numpy.array([1E-7, 2E7])
    xs = numpy.array([11.4])
    u235.setMultigroupElasticXS(xs_energies, xs)

    # Set one group potential elastic scattering xs for u-238
    xs = numpy.array([11.3])
    u238.setMultigroupElasticXS(xs_energies, xs)

    py_printf('INFO', 'Initializing fuel and moderator materials...')
    
    # Define moderator material
    moderator = Material('moderator')
    moderator.setDensity(0.7, 'g/cc')
    moderator.addIsotope(h1, 2.0)
    moderator.addIsotope(o16, 1.0)

    # Define fuel material
    fuel = Material('fuel')
    fuel.setDensity(10.2, 'g/cc')
    fuel.addIsotope(u235, 0.03035)
    fuel.addIsotope(u238, 0.9695)
    fuel.addIsotope(o16, 2.0)
    
    py_printf('INFO', 'Initializing fuel and moderator regions...')
    
    # Define moderator region
    region_mod = Region('moderator', MODERATOR)
    region_mod.setMaterial(moderator)
    region_mod.setFuelRadius(0.4096)
    region_mod.setPitch(1.26)
        
    # Define fuel region
    region_fuel = Region('fuel', FUEL)
    region_fuel.setMaterial(fuel)
    region_fuel.setFuelRadius(0.4096)
    region_fuel.setPitch(1.26)
    
    py_printf('INFO', 'Initializing the geometry...')

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(HOMOGENEOUS_EQUIVALENCE)
    geometry.addRegion(region_mod)
    geometry.addRegion(region_fuel)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)
    geometry.setDancoffFactor(dancoff)

    py_printf('INFO', 'Initializing flux tallies...')

    # Create Tallies for the fluxes
    total_flux = TallyFactory.createTally(geometry, FLUX, 'total')
    moderator_flux = TallyFactory.createTally(moderator, FLUX, 'moderator')
    fuel_flux = TallyFactory.createTally(fuel, FLUX, 'fuel')
    total_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
    moderator_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
    fuel_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)

	# Register tallies
    TallyBank.registerTally(total_flux)
    TallyBank.registerTally(moderator_flux)
    TallyBank.registerTally(fuel_flux)

	# Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

    py_printf('INFO', 'Writing tally batch statistics to output file...')

	# Dump batch statistics to output files to some new directory
    TallyBank.outputBatchStatistics()

    py_printf('INFO', 'Plotting fluxes...')

    # Plotting xs, flux, thermal scattering
    plotter.plotFluxes([total_flux, moderator_flux, fuel_flux])

    py_printf('INFO', 'Plotting microscopic and macroscopic cross-sections...')

    plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission'])
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    plotter.plotMacroXS(fuel, ['capture', 'elastic', 'fission'])
    plotter.plotMacroXS(moderator, ['capture', 'elastic', 'fission'])

    py_printf('TITLE', 'Finished')


if __name__ == '__main__':
    
    main()

