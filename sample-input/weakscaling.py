import time
import matplotlib.pyplot as plt
import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
from pinspec.log import *

def main():

    # Set main simulation params
    num_neutrons_per_batch = 100000
    setOutputDirectory('weak-scaling')

    log_setlevel(INFO)

    py_printf('TITLE', 'Starting a weak scaling multi-threading study')

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
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)

    py_printf('INFO', 'Initializing flux tally...')

    # Create a tally for the flux
    flux = TallyFactory.createTally(region_mix, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC) 
    
	# Register the tallies
    TallyBank.registerTally(flux)

    runtimes = numpy.zeros(12)
    threads = numpy.linspace(1,12,12)

    # Run simulation with 1-13 threads
    for num_threads in threads: 

        geometry.setNumThreads(int(num_threads))
        geometry.setNumBatches(2*int(num_threads))

	    # Run Monte Carlo simulation
        start_time = time.time()
        geometry.runMonteCarloSimulation();
        end_time = time.time()
        runtimes[num_threads-1] = (end_time - start_time)

    # Plot the runtime vs. thread count
    fig = plt.figure()
    plt.plot(threads, runtimes)
    plt.title('Weak Scaling: Runtime vs. Thread Count')
    plt.xlabel('# Threads')
    plt.ylabel('Runtime [sec]')
    plt.grid()
    plt.savefig(getOutputDirectory() + '/threadscaling.png')


if __name__ == '__main__':

    main()  

