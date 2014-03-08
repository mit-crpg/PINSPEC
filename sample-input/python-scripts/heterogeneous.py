import numpy
import math
import random
import matplotlib.pyplot as plt
from pinspec import *
import pinspec.process as process
import pinspec.plotter as plotter
from pinspec.log import *


###############################################################################
##    This python script uses PINSPEC to simulate the flux in fuel and 
##    moderator regions using ray tracing across a heterogeneous geometry.
###############################################################################


###############################################################################
###########################  Main Simulation Paramters  #######################
###############################################################################

# Simulation parameters
num_batches = 25
num_neutrons_per_batch = 10000
num_threads = 4
nu = 2.45
set_output_directory('heterogeneous')

# Geometric parameters
radius_fuel = 0.4096
radius_gap = 0.4178
radius_clad = 0.4750
pitch = 1.26
pitch_squared = pitch**2.0
dancoff = 0.277
source_sampling_radius = 2.0      # sphere within which to do rejection sampling

py_set_log_level('INFO')

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
###############################  Create Surfaces ##############################
###############################################################################

py_printf('INFO', 'Initializing bounding surfaces...')
left = XPlane('left edge')
right = XPlane('right edge')
bottom = YPlane('bottom edge')
top = YPlane('top edge')
pin = ZCylinder('fuel boundary')

left.setX(-pitch/2.0)
right.setX(pitch/2.0)
bottom.setY(-pitch/2.0)
top.setY(pitch/2.0)
pin.setRadius(radius_fuel)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)
pin.setBoundaryType(INTERFACE)


###############################################################################
###############################  Create Regions ###############################
###############################################################################

py_printf('INFO', 'Initializing fuel and moderator regions...')
region_mod = BoundedModeratorRegion('moderator')
region_mod.setMaterial(moderator)
region_mod.addBoundingSurface(+1, left)
region_mod.addBoundingSurface(-1, right)
region_mod.addBoundingSurface(+1, bottom)
region_mod.addBoundingSurface(-1, top)
region_mod.addBoundingSurface(+1, pin)
region_mod.setVolume(vol_moderator)

region_fuel = BoundedFuelRegion('fuel')
region_fuel.setMaterial(fuel)
region_fuel.addBoundingSurface(-1, pin)
region_fuel.setVolume(vol_fuel)

###############################################################################
###############################  Create Geometry ##############################
###############################################################################

py_printf('INFO', 'Initializing the geometry...')
geometry = Geometry(HETEROGENEOUS)
geometry.addRegion(region_mod)
geometry.addRegion(region_fuel)
geometry.setNumBatches(num_batches)
geometry.setNeutronsPerBatch(num_neutrons_per_batch)
geometry.setNumThreads(num_threads)
geometry.setSourceSamplingRadius(0.5)


###############################################################################
################################  Create Tallies ##############################
###############################################################################

py_printf('INFO', 'Initializing tallies...')

# Create Tallies for the energy high-resolution fluxes for plotting
total_flux = TallyFactory.createTally(geometry, FLUX, 'total')
moderator_flux = TallyFactory.createTally(region_mod, FLUX, 'moderator')
fuel_flux = TallyFactory.createTally(region_fuel, FLUX, 'fuel')
total_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
moderator_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
fuel_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)

# Create tallies for computing k-inf
tot_fiss_rate = TallyFactory.createTally(geometry, FISSION_RATE)
tot_abs_rate = TallyFactory.createTally(geometry, ABSORPTION_RATE)
tot_fiss_rate.generateBinEdges(0.0, 1E7, 1, EQUAL)
tot_abs_rate.generateBinEdges(0.0, 1E7, 1, EQUAL)

# Register tallies
TallyBank.registerTally(total_flux)
TallyBank.registerTally(moderator_flux)
TallyBank.registerTally(fuel_flux)
TallyBank.registerTally(tot_fiss_rate)
TallyBank.registerTally(tot_abs_rate)


###############################################################################
#########################  Run Monte Carlo Simulation #########################
###############################################################################

# Run Monte Carlo simulation
geometry.runMonteCarloSimulation()


###############################################################################
##################################  Plotting  #################################
###############################################################################

py_printf('INFO', 'Plotting the geometry...')

# Plot xy slices
plotter.plotSlice(geometry, plane='XY', lim1=[-1., 1.], lim2=[-1., 1.], \
                    filename='fuel', gridsize=50)
plotter.plotSlice(geometry, plane='XY', lim1=[-1., 1.], lim2=[-1., 1.], \
                    filename='moderator', gridsize=50)

# Plot xz slices
plotter.plotSlice(geometry, plane='XZ', lim1=[-1., 1.], lim2=[-1., 1.], \
                    filename='fuel', gridsize=50)
plotter.plotSlice(geometry, plane='XZ', lim1=[-1., 1.], lim2=[-1., 1.], \
                    filename='moderator', gridsize=50)

# Plot yz slices
plotter.plotSlice(geometry, plane='YZ', lim1=[-1., 1.], lim2=[-1., 1.], \
                    filename='fuel', gridsize=50)
plotter.plotSlice(geometry, plane='YZ', lim1=[-1., 1.], lim2=[-1., 1.], \
                    filename='moderator', gridsize=50)

py_printf('INFO', 'Plotting the fission site source distribution...')
plotter.plotFissionSourceDist(geometry)


py_printf('INFO', 'Plotting neutron tracks...')

plotter.trackANeutron(geometry, num_moves=100, plane='xy', \
                          lim1=[-0.75, 0.75], lim2=[-0.75,0.75])
plotter.trackANeutron(geometry, num_moves=100, plane='yz', \
                          lim1=[-0.75, 0.75], lim2=[-3., 3.])
plotter.trackANeutron(geometry, num_moves=100, plane='xz', \
                          lim1=[-0.75, 0.75], lim2=[-3.,3.])

# Plot the fluxes
plotter.plotFlux([total_flux, moderator_flux, fuel_flux], uselegend=True)


# Compute k-infinity using tally arithmetic operators and print to screen
k_inf = tot_fiss_rate * nu  / tot_abs_rate
k_inf.setTallyName('k-infinity')
k_inf.printTallies(uncertainties=True)

py_printf('TITLE', 'Finished')
