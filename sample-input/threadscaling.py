import time
import matplotlib.pyplot as plt
import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
from pinspec.log import *

def main():

    # Set main simulation params
    num_batches = 25
    num_neutrons_per_batch = 200000
    setlevel(INFO)

    # Call SLBW to create XS
    filename = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
    T=300 #Temp in Kelvin of target nucleus
    SLBW.SLBWXS(filename,T)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    zr90 = Isotope('Zr-90')
 
    # Define materials
    mix = Material('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(u238, 0.40)
    mix.addIsotope(u235, .02)
    mix.addIsotope(zr90, 0.16)    

    runtimes = numpy.zeros(12)
    threads = numpy.linspace(1,12,12)

    # Run simulation with 1-13 threads
    for num_threads in threads: 
           
        # Define regions
        region_mix = Region('infinite medium', INFINITE)
        region_mix.setMaterial(mix)

        # Create a tally for the flux
        flux = TallyFactory.createTally('total flux', region_mix, FLUX)
        flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC)

        abs_rate = TallyFactory.createTally('absorption rate', region_mix, ABSORPTION_RATE)
        abs_rate_bin_edges = numpy.array([0.1, 1., 5., 10., 100., 1000.])
        abs_rate.setBinEdges(abs_rate_bin_edges)

		TallyBank.registerTally(abs_rate)
		TallyBank.registerTally(flux)

        # Define geometry
        geometry = Geometry()
        geometry.setSpatialType(INFINITE_HOMOGENEOUS)
        geometry.addRegion(region_mix)
        geometry.setNumBatches(num_batches)
        geometry.setNeutronsPerBatch(num_neutrons_per_batch)
        geometry.setNumThreads(int(num_threads))

	    # Run Monte Carlo simulation
        start_time = time.time()
        geometry.runMonteCarloSimulation();
        end_time = time.time()
        runtimes[num_threads-1] = (end_time - start_time)


    # Plot the runtime vs. thread count
    fig = plt.figure()
    print 'num_threads = ' + str(threads.size) + ' runtimes = ' + str(runtimes.size)
    plt.plot(threads, runtimes)
    print 'threads = ' + str(threads)
    print 'runtimes = ' + str(runtimes)
    plt.title('Runtime vs. Thread Count')
    plt.xlabel('# Threads')
    plt.ylabel('Runtime [sec]')
    plt.grid()
    plt.savefig('threadscaling.png')


if __name__ == '__main__':

    main()  

