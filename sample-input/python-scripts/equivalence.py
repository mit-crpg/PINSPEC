import numpy
import math
from pinspec import *
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

###############################################################################
##    This python script uses PINSPEC to generate results for homework 4 of
##    the 22.211 Introduction to Reactor Physics course taught at MIT by
##    Prof. Kord Smith in spring 2012.
###############################################################################


###############################################################################
###########################  Main Simulation Paramters  #######################
###############################################################################

# Simulation parameters
num_batches = 50
num_neutrons_per_batch = 10000
num_threads = 4
nu = 2.45
setOutputDirectory('equivalence')

# Geometric parameters
radius_fuel = 0.4096
radius_gap = 0.4178
radius_clad = 0.4750
pitch = 1.26
pitch_squared = pitch**2.0
dancoff = 0.277

py_setlevel('INFO')

py_printf('TITLE', 'Starting a heterogeneous-homogeneous equivalence ' + \
              'calculation')


###############################################################################
###################  Compute Homogenized Number Densities  ####################
###############################################################################

py_printf('INFO', 'Computing homogenized number densities...')

# 2-region homogenized densities [g/cc] and enrichment
rho_fuel = 10.2
rho_clad = 6.54
rho_coolant = 0.9966
enrichment = 0.03035

# Compute material number densities
N_A = 6.023E23
N_u238 = rho_fuel * N_A * (1.0 - enrichment) / ((238.0 *
                    (1.0 - enrichment)) + (235.0*enrichment) + (16.0*2.0))
N_u235 = rho_fuel * N_A * enrichment / ((238.0 *
		    (1.0 - enrichment)) + (235.0*enrichment) + (16.0*2.0))
N_o16 = rho_fuel * N_A * 2.0 / ((238.0 *
		    (1.0 - enrichment)) + (235.0*enrichment) + (16.0*2.0))
N_zr90 = rho_clad * N_A / 90.0
N_h2o = rho_coolant * N_A / 18.0
N_h1 = rho_coolant * N_A * 2.0 / 18.0

# Compute 2-region pin cell volumes [cm^3]
vol_fuel = math.pi * radius_fuel**2.0
vol_gap = math.pi * (radius_gap**2.0 - radius_fuel**2.0)
vol_clad = math.pi * (radius_clad**2.0 - radius_gap**2.0)
vol_coolant = pitch_squared	 - math.pi * radius_clad**2.0
vol_moderator = vol_gap + vol_clad + vol_coolant;
vol_total = vol_fuel + vol_moderator;

#  Compute homogenized moderator number densities using volume weighting
N_h2o *= (vol_coolant / vol_moderator);
N_h1 *= (vol_coolant / vol_moderator);
N_zr90 *= (vol_clad / vol_moderator);
N_moderator = N_h2o + N_zr90 + N_h1
N_fuel = N_u238 + N_u235 + N_o16

rho_moderator = (vol_coolant / vol_moderator) * rho_coolant + \
                  (vol_clad / vol_moderator) * rho_clad


###############################################################################
##############################  Create Isotopes ###############################
###############################################################################

py_printf('INFO', 'Initializing isotopes...')
h1 = Isotope('H-1')
o16 = Isotope('O-16')
u235 = Isotope('U-235')
u238 = Isotope('U-238')
zr90 = Isotope('Zr-90')

# Neglect thermal scattering in O-16, U-235, U-238
o16.neglectThermalScattering()
u235.neglectThermalScattering()
u238.neglectThermalScattering()
zr90.neglectThermalScattering()

# Set one group potential elastic scattering xs for u-235
xs_energies = numpy.array([1E-7, 2E7])
xs = numpy.array([11.4])
u235.setMultigroupElasticXS(xs_energies, xs)

# Set one group potential elastic scattering xs for u-238
xs = numpy.array([11.3])
u238.setMultigroupElasticXS(xs_energies, xs)

# Zero out capture for Zr-90 and O-16
xs = numpy.array([0.0])
zr90.setMultigroupCaptureXS(xs_energies, xs)
o16.setMultigroupCaptureXS(xs_energies, xs)


###############################################################################
##############################  Create Materials ##############################
###############################################################################

py_printf('INFO', 'Initializing fuel and moderator materials...')
moderator = Material('moderator')
moderator.setDensity(rho_moderator, 'g/cc')
moderator.addIsotope(h1, N_h1 / N_moderator)
moderator.addIsotope(o16, N_h2o / N_moderator)
moderator.addIsotope(zr90, N_zr90 / N_moderator)

# Define fuel material
fuel = Material('fuel')
fuel.setDensity(rho_fuel, 'g/cc')
fuel.addIsotope(u235, N_u235 / N_fuel)
fuel.addIsotope(u238, N_u238 / N_fuel)
fuel.addIsotope(o16, N_o16 / N_fuel)


###############################################################################
###############################  Create Regions ###############################
###############################################################################

py_printf('INFO', 'Initializing fuel and moderator regions...')
region_mod = EquivalenceModeratorRegion('moderator')
region_mod.setMaterial(moderator)
        
region_fuel = EquivalenceFuelRegion('fuel')
region_fuel.setMaterial(fuel)


###############################################################################
###############################  Create Geometry ##############################
###############################################################################
    
py_printf('INFO', 'Initializing the geometry...')
geometry = Geometry(HOMOGENEOUS_EQUIVALENCE)
geometry.addRegion(region_mod)
geometry.addRegion(region_fuel)
geometry.setFuelPinRadius(radius_fuel)
geometry.setPinCellPitch(pitch)
geometry.setNumBatches(num_batches)
geometry.setNeutronsPerBatch(num_neutrons_per_batch)
geometry.setNumThreads(num_threads)
geometry.setDancoffFactor(dancoff)


###############################################################################
################################  Create Tallies ##############################
###############################################################################

py_printf('INFO', 'Initializing tallies...')

# Create Tallies for the energy high-resolution fluxes for plotting
total_flux = TallyFactory.createTally(geometry, FLUX, 'total')
moderator_flux = TallyFactory.createTally(moderator, FLUX, 'moderator')
fuel_flux = TallyFactory.createTally(fuel, FLUX, 'fuel')
total_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
moderator_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
fuel_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)

# Create Tallies for moderator-to-fuel flux ratios
fuel_flux_ratio = TallyFactory.createTally(fuel, FLUX, 'Fuel Flux')
moderator_flux_ratio = TallyFactory.createTally(moderator,FLUX,'Moderator Flux')
flux_bin_edges = numpy.array([0., 0.1, 0.5, 1., 6., 10., 25., \
                              50., 1E2, 1E3, 1E4, 1E5, 5E5, 1E7])
moderator_flux_ratio.setBinEdges(flux_bin_edges)
fuel_flux_ratio.setBinEdges(flux_bin_edges)

# Create tallies for two group cross-sections
total_flux_xs = TallyFactory.createTally(geometry, FLUX)
elastic_rate = TallyFactory.createTally(geometry, ELASTIC_RATE)
capture_rate = TallyFactory.createTally(geometry, CAPTURE_RATE)
fission_rate = TallyFactory.createTally(geometry, FISSION_RATE)
absorb_rate = TallyFactory.createTally(geometry, ABSORPTION_RATE)
transport_rate = TallyFactory.createTally(geometry, TRANSPORT_RATE)
diffusion_rate = TallyFactory.createTally(geometry, DIFFUSION_RATE)
total_rate = TallyFactory.createTally(geometry, COLLISION_RATE)

group_xs_edges = numpy.array([0., 0.625, 1E7])
total_flux_xs.setBinEdges(group_xs_edges)
capture_rate.setBinEdges(group_xs_edges)
fission_rate.setBinEdges(group_xs_edges)
absorb_rate.setBinEdges(group_xs_edges)
elastic_rate.setBinEdges(group_xs_edges)
transport_rate.setBinEdges(group_xs_edges)
diffusion_rate.setBinEdges(group_xs_edges)
total_rate.setBinEdges(group_xs_edges)

# Create tallies for computing k-inf
tot_fiss_rate = TallyFactory.createTally(geometry, FISSION_RATE)
tot_abs_rate = TallyFactory.createTally(geometry, ABSORPTION_RATE)
tot_fiss_rate.generateBinEdges(0.0, 1E7, 1, EQUAL)
tot_abs_rate.generateBinEdges(0.0, 1E7, 1, EQUAL)

# Register tallies
TallyBank.registerTally(total_flux)
TallyBank.registerTally(moderator_flux)
TallyBank.registerTally(fuel_flux)
TallyBank.registerTally(fuel_flux_ratio)
TallyBank.registerTally(moderator_flux_ratio)
TallyBank.registerTally(total_flux_xs)
TallyBank.registerTally(elastic_rate)
TallyBank.registerTally(capture_rate)
TallyBank.registerTally(fission_rate)
TallyBank.registerTally(absorb_rate)
TallyBank.registerTally(transport_rate)
TallyBank.registerTally(diffusion_rate)
TallyBank.registerTally(total_rate)
TallyBank.registerTally(tot_fiss_rate)
TallyBank.registerTally(tot_abs_rate)


###############################################################################
#########################  Run Monte Carlo Simulation #########################
###############################################################################

# Run Monte Carlo simulation
geometry.runMonteCarloSimulation();


###############################################################################
############################  Process Output Data #############################
###############################################################################

py_printf('INFO', 'Plotting microscopic and macroscopic cross-sections...')
plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission'])
plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'])
plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
plotter.plotMacroXS(fuel, ['capture', 'elastic', 'fission'])
plotter.plotMacroXS(moderator, ['capture', 'elastic', 'fission'])

py_printf('INFO', 'Computing group cross-sections...')
elastic_xs = process.GroupXS(total_flux_xs, elastic_rate)
capture_xs = process.GroupXS(total_flux_xs, capture_rate)
fission_xs = process.GroupXS(total_flux_xs, fission_rate)
absorb_xs = process.GroupXS(total_flux_xs, absorb_rate)
transport_xs = process.GroupXS(total_flux_xs, transport_rate)
diffusion_coeff = process.GroupXS(total_flux_xs, diffusion_rate)
total_xs = process.GroupXS(total_flux_xs, total_rate)

elastic_xs.printXS(uncertainties=True)
capture_xs.printXS(uncertainties=True)
fission_xs.printXS(uncertainties=True)
absorb_xs.printXS(uncertainties=True)
transport_xs.printXS(uncertainties=True)
diffusion_coeff.printXS(uncertainties=True)
total_xs.printXS(uncertainties=True)

# Compute k-infinity using tally arithmetic operators and print to screen
k_inf = tot_fiss_rate * nu  / tot_abs_rate
k_inf.setTallyName('k-infinity')
k_inf.printTallies(uncertainties=True)

# Compute moderator-to-fuel flux ratios and print to screen
flux_ratio = moderator_flux_ratio / fuel_flux_ratio
flux_ratio.setTallyName('Moderator-to-Fuel Flux Ratios')
flux_ratio.printTallies(uncertainties=True)

# Plot the fluxes
plotter.plotFlux([total_flux, moderator_flux, fuel_flux])

py_printf('TITLE', 'Finished')
