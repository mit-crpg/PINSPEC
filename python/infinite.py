import matplotlib.pyplot as plt
import numpy as np
from pinspec import *
from numpy import *
#from plotter import *


def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 10000
    num_threads = 4
    log_setlevel(INFO)


    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')


    ############################################################################
    # EXAMPLE: How to retrieve micro xs from C++ and plot them
    ############################################################################
    # First, retrieve the xs and energies
    u235_num_xs = u235.getNumXSEnergies()
    u235_capture_xs = u235.retrieveXS(u235_num_xs, 'capture')
    u235_elastic_xs = u235.retrieveXS(u235_num_xs, 'elastic')
    u235_fission_xs = u235.retrieveXS(u235_num_xs, 'fission')
    u235_absorption_xs = u235.retrieveXS(u235_num_xs, 'absorption')
    u235_total_xs = u235.retrieveXS(u235_num_xs, 'total')
    #Note: after rescaling, all xs types have the same energy grid so only
    #need to retrieve it for one xs type
    u235_xs_energies = u235.retrieveXSEnergies(u235_num_xs)    

    # Plot all of the xs on the same scale
    fig = plt.figure()
    plt.plot(u235_xs_energies, u235_capture_xs, lw=1)
    plt.plot(u235_xs_energies, u235_elastic_xs, lw=1)
    plt.plot(u235_xs_energies, u235_fission_xs, lw=1)
    plt.plot(u235_xs_energies, u235_absorption_xs, lw=1)
    plt.plot(u235_xs_energies, u235_total_xs, lw=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [ev]')
    plt.ylabel('Micro XS [barns]')
    plt.title(u235.getIsotopeType() + ' Micro XS')
    plt.legend(['capture', 'elastic', 'fission', 'absorption', 'total'])
    fig.savefig(u235.getIsotopeType() + '_micro_xs.png')

    
    # Define materials
    mix = Material()
    mix.setMaterialName('fuel moderator mix')
    mix.setDensity(10., 'g/cc')
    mix.addIsotope(h1, 0.1)
    mix.addIsotope(o16, 0.05)
    mix.addIsotope(u235, 0.03)
    mix.addIsotope(u238, 0.97)

    log_printf(INFO, "Added isotopes")


    ############################################################################
    #EXAMPLE: How to retrieve macro xs from C++ and plot them
    ############################################################################
    # First, retrieve the xs and energies 
    mix_num_xs = mix.getNumXSEnergies()
    mix_capture_xs = mix.retrieveXS(mix_num_xs, 'capture')
    mix_elastic_xs = mix.retrieveXS(mix_num_xs, 'elastic')
    mix_fission_xs = mix.retrieveXS(u235_num_xs, 'fission')
    mix_absorption_xs = mix.retrieveXS(u235_num_xs, 'absorption')
    mix_total_xs = mix.retrieveXS(u235_num_xs, 'total')
    #Note: after rescaling, all xs types have the same energy grid so only
    #need to retrieve it for one xs type
    mix_xs_energies = mix.retrieveXSEnergies(mix_num_xs)    

    # Plot all of the xs on the same scale
    fig = plt.figure()
    plt.plot(mix_xs_energies, mix_capture_xs, lw=1)
    plt.plot(mix_xs_energies, mix_elastic_xs, lw=1)
    plt.plot(mix_xs_energies, mix_fission_xs, lw=1)
    plt.plot(mix_xs_energies, mix_absorption_xs, lw=1)
    plt.plot(mix_xs_energies, mix_total_xs, lw=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [ev]')
    plt.ylabel('Macro XS [cm^-1]')
    plt.title(mix.getMaterialName() + ' Macro XS')
    plt.legend(['capture', 'elastic', 'fission', 'absorption', 'total'])
    fig.savefig(mix.getMaterialName() + '_macro_xs.png')


    
    # Define regions
    region_mix = Region()
    region_mix.setRegionName('infinite medium fuel/moderator mix')
    region_mix.setRegionType(INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, "Made mixture region")
    

	# Define tallies - give them to Regions, Materials, or Isotopes
	# This part is really where we need to know how to pass float
    # arrays to/from SWIG

    # Create a tally for the flux
    flux = Tally('total flux', REGION, FLUX)
    flux.generateBinEdges(1E-7, 1E7, 10000, LOGARITHMIC)
    region_mix.addTally(flux)

    ############################################################################
    #EXAMPLE: How to set tally bin edges 
    ############################################################################
    # Create a tally for the absorption rate
    abs_rate = Tally('absorption rate', MATERIAL, ABSORPTION_RATE)
    abs_rate_bin_edges = np.array([0.1, 1., 5., 10., 100., 1000.])
    abs_rate.setBinEdges(abs_rate_bin_edges)
    mix.addTally(abs_rate)

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    log_printf(INFO, "Made geometry")

	# Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

        
    ############################################################################
    #EXAMPLE: How to retrieve tally data. 
    ############################################################################
    # Flux
    num_bins = flux.getNumBins()
    flux.computeBatchStatistics()
    flux_bin_centers = flux.retrieveTallyCenters(num_bins)
    flux_mu = flux.retrieveTallyMu(num_bins)
    flux_variance = flux.retrieveTallyVariance(num_bins)
    flux_std_dev = flux.retrieveTallyStdDev(num_bins)
    flux_rel_err = flux.retrieveTallyRelErr(num_bins)

    # Plot the flux
    fig = plt.figure()
    plt.plot(flux_bin_centers, flux_mu, lw=1)
    plt.xscale('log')
    plt.xlabel('Energy [ev]')
    if (flux.getTallyType() == FLUX):
        plt.ylabel('Flux')
    plt.title(flux.getTallyName() + ' average')
    fig.savefig(flux.getTallyName() + '_average.png')

    # Absorption rate step function
    num_bins = abs_rate.getNumBins()
    abs_rate.computeBatchStatistics()
    abs_rate_bin_centers = abs_rate.retrieveTallyCenters(num_bins)
    abs_rate_mu = abs_rate.retrieveTallyMu(num_bins)
    abs_rate_variance = abs_rate.retrieveTallyVariance(num_bins)
    abs_rate_std_dev = abs_rate.retrieveTallyStdDev(num_bins)
    abs_rate_rel_err = abs_rate.retrieveTallyRelErr(num_bins)

    # Plot the absorption rate
    fig = plt.figure()
    plt.step(abs_rate_bin_centers, abs_rate_mu, lw=1)
    plt.xscale('log')
    plt.xlabel('Energy [ev]')
    if (abs_rate.getTallyType() == ABSORPTION_RATE):
        plt.ylabel('Absorption Rate')
    plt.title(abs_rate.getTallyName() + ' average')
    fig.savefig(abs_rate.getTallyName() + '_average.png')

	# Dump batch statistics to output files to some new directory - gives segmentation fault right now
    #geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')


if __name__ == '__main__':
    
    main()
