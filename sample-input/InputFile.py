from pinspec import *
from pinspec.log import *
import pinspec.slbw as slbw
import pinspec.plotter as plotter
import pinspec.process as process
import numpy

py_setlevel('NORMAL')
py_printf('TITLE', 'Simulating an infinite medium homogenized problem')

# Call SLBW to broaden XS
py_printf('NORMAL', 'Creating SLBW xs')

slbw.SLBWXS('U-238',1200,'capture')
slbw.generatePotentialScattering('U-238')
slbw.compareXS('U-238', XStype='capture')

# Define isotopes
py_printf('NORMAL', 'Initializing isotopes...')

u235 = Isotope('U-235')
u238 = Isotope('U-238')
o16 = Isotope('O-16')
h1 = Isotope('H-1')

# Define the material
py_printf('NORMAL', 'Initializing homogeneous fuel-moderator mix material...')

mix = Material('Fuel Moderator Homogeneous Mix')
mix.setDensity(5., 'g/cc')
mix.addIsotope(o16, 1.0)
mix.addIsotope(h1, 1.0)
mix.addIsotope(u238, 0.01)
mix.addIsotope(u235, .0025)

#Define the region
py_printf('NORMAL', 'Initializing infnite homogeneous medium region...')

region_mix = InfiniteMediumRegion('Infinite Medium')
region_mix.setMaterial(mix)

#Define the geometry
py_printf('NORMAL', 'Initializing the geometry...')

geometry = Geometry(INFINITE_HOMOGENEOUS)
geometry.addRegion(region_mix)
num_batches = 10
geometry.setNumBatches(num_batches)
num_neutrons_per_batch = 10000
geometry.setNeutronsPerBatch(num_neutrons_per_batch)
num_threads = 12
geometry.setNumThreads(num_threads)
geometry.setSourceSamplingRadius(0.5)
setOutputDirectory('infinite_homogeneous_output')

#Define the tallies
py_printf('NORMAL', 'Initializing flux tally...')

flux = TallyFactory.createTally(region_mix, FLUX, 'Flux in Infinite Region')
elastic_rate = TallyFactory.createTally(geometry, ELASTIC_RATE)
capture_rate = TallyFactory.createTally(geometry, CAPTURE_RATE)
fission_rate = TallyFactory.createTally(geometry, FISSION_RATE)
TallyBank.registerTally(flux)
TallyBank.registerTally(elastic_rate)
TallyBank.registerTally(capture_rate)
TallyBank.registerTally(fission_rate)

#Run MC simulation
geometry.runMonteCarloSimulation()
TallyBank.outputBatchStatistics()

#plots
plotter.plotFlux(flux, loglog=True, uselegend=False, filename='flux', title='Infinite Homogeneous Medium flux')
plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'])
plotter.plotThermalScattering(h1)
