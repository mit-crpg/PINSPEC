import numpy
import time
import matplotlib.pyplot as plt
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter


def main():

    # Set main simulation params
    num_batches = 25
    num_neutrons_per_batch = 200000
    log_setlevel(INFO)

    setXSLibDirectory('../xs-lib/')   # This is also a default, but set it as example

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
    mix = Material()
    mix.setMaterialName('Fuel Moderator Mix')
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
        flux = Tally('total flux', GEOMETRY, FLUX)
        flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC)

        ############################################################################
        #EXAMPLE: How to set tally bin edges 
        ############################################################################
        # Create a tally for the absorption rate
        abs_rate = Tally('absorption rate', REGION, ABSORPTION_RATE)
        abs_rate_bin_edges = numpy.array([0.1, 1., 5., 10., 100., 1000.])
        abs_rate.setBinEdges(abs_rate_bin_edges)
        region_mix.addTally(abs_rate)

        # Define geometry
        geometry = Geometry()
        geometry.setSpatialType(INFINITE_HOMOGENEOUS)
        geometry.addRegion(region_mix)
        geometry.setNumBatches(num_batches)
        geometry.setNeutronsPerBatch(num_neutrons_per_batch)
        geometry.setNumThreads(int(num_threads))
        geometry.addTally(flux)

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

