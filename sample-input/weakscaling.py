import time
import matplotlib.pyplot as plt
import numpy
from pinspec import *
import pinspec.plotter as plotter
from pinspec.log import *


###############################################################################
##    This python script uses PINSPEC to generate results for a weak scaling
##    study of PINSPEC's parallel performance. This file runs a series 
##    of infinite medium spectral calculations and varies the number of
##    shared memory (OpenMP) parallel threads with the number of batches.
###############################################################################


###############################################################################
###########################  Main Simulation Paramters  #######################
###############################################################################

num_neutrons_per_batch = 10000
setOutputDirectory('weakscaling')
py_setlevel('INFO')

py_printf('TITLE', 'Starting a weak scaling multi-threading study')


###############################################################################
##############################  Create Isotopes ###############################
###############################################################################

py_printf('INFO', 'Initializing isotopes...')
h1 = Isotope('H-1')
b10 = Isotope('B-10')
o16 = Isotope('O-16')
u235 = Isotope('U-235')
u238 = Isotope('U-238')
zr90 = Isotope('Zr-90')    


###############################################################################
##############################  Create Materials ##############################
###############################################################################

py_printf('INFO', 'Initializing fuel-moderator mix material...')
mix = Material('Fuel Moderator Mix')
mix.setDensity(5., 'g/cc')
mix.addIsotope(b10, .0000001)
mix.addIsotope(o16, 1.0)
mix.addIsotope(h1, 1.0)
mix.addIsotope(u238, 0.01)
mix.addIsotope(u235, .0025)
mix.addIsotope(zr90, .16)


###############################################################################
###############################  Create Regions ###############################
###############################################################################

py_printf('INFO', 'Initializing fuel-moderator mix region...')
region_mix = InfiniteMediumRegion('infinite medium')
region_mix.setMaterial(mix)


###############################################################################
###############################  Create Geometry ##############################
###############################################################################
    
py_printf('INFO', 'Initializing the geometry...')
geometry = Geometry(INFINITE_HOMOGENEOUS)
geometry.addRegion(region_mix)
geometry.setNeutronsPerBatch(num_neutrons_per_batch)


###############################################################################
################################  Create Tallies ##############################
###############################################################################

py_printf('INFO', 'Initializing flux tally...')
flux = TallyFactory.createTally(region_mix, FLUX)
flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC) 

# Register the tallies
TallyBank.registerTally(flux)


###############################################################################
#########################  Run Monte Carlo Simulation #########################
###############################################################################

runtimes = numpy.zeros(12)
threads = numpy.linspace(1,12,12)

# Run simulation with 1-12 threads and 2-24 batches
for num_threads in threads: 

    geometry.setNumThreads(int(num_threads))
    geometry.setNumBatches(2*int(num_threads))

    # Run Monte Carlo simulation
    start_time = time.time()
    geometry.runMonteCarloSimulation();
    end_time = time.time()
    runtimes[num_threads-1] = (end_time - start_time)


###############################################################################
############################  Process Output Data #############################
###############################################################################

# Plot the runtime vs. thread count
fig = plt.figure()
plt.plot(threads, runtimes)
plt.title('Weak Scaling: Runtime vs. Thread Count')
plt.xlabel('# Threads')
plt.ylabel('Runtime [sec]')
plt.grid()
plt.savefig(getOutputDirectory() + '/threadscaling.png')

py_printf('HEADER', 'Finished')
