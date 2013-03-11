import matplotlib.pyplot as plt
import numpy as np
from pinspec import *

# Series of functions use to plot XS data (energy, xs) that is saved
# in text files output by pinspec. The plotter creates
# a directory to store the text files, reads the text files
# and outputs png plots, and deletes the directory
# used to store the text files.

# Function to plot the microscopic cross section for a
# given isotope and array of reactions
def plotMicroXS(isotope, rxns):
    
    # set input and output file names
    filename = isotope.getIsotopeType() + '-micro-xs.png'
    
	# retrieve xs energies
    num_energies = isotope.getNumXSEnergies()
    energies = isotope.retrieveXSEnergies(num_energies)

    # make figure
    fig = plt.figure()
        
    # loop over rxns and plot
    for rxn in rxns:
        xs = isotope.retrieveXS(num_energies, rxn)
        plt.plot(energies, xs, lw=1)
        
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [ev]')
    plt.ylabel('Micro XS [barns]')
    plt.title(isotope.getIsotopeType() + ' Micro XS')
    plt.legend(rxns)
    fig.savefig(filename)

# function to plot the macroscopic cross section for a
# given material and arrary of reactions
def plotMacroXS(material, rxns):
    
    # set input and output file names
    filename = material.getMaterialName() + '-macro-xs.png'
    
	# retrieve xs energies
    num_energies = material.getNumXSEnergies()
    energies = material.retrieveXSEnergies(num_energies)
    
    # make figure
    fig = plt.figure()
    
    # loop over rxns and plot
    for rxn in rxns:
        xs = material.retrieveXS(num_energies, rxn)
        plt.plot(energies, xs, lw=1)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [ev]')
    plt.ylabel('Macro XS [cm^-1]')
    plt.title(material.getMaterialName() + ' Macro XS')
    plt.legend(rxns)
    fig.savefig(filename)

# Function that plots the flux spectrum
def plotFlux(flux):

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
    plt.yscale('log')
    plt.xlabel('Energy [ev]')
    if (flux.getTallyType() == FLUX):
        plt.ylabel('Flux')

    plt.title(flux.getTallyName() + ' average')
    fig.savefig(flux.getTallyName() + '_average.png')


# Function to plot the fission CDF and sample from the
# fission CDF to generate a fission spectrum
def plotFissionSpectrum():
    
    # Generate fission CDF
    fissioner = Fissioner()
    fissioner.setNumBins(10000)
    fissioner.setEMax(20)
    fissioner.buildCDF()
    cdf = fissioner.retrieveCDF(fissioner.getNumBins())
    cdf_energies = fissioner.retrieveCDFEnergies(fissioner.getNumBins())

    # Plot fission CDF
    fig = plt.figure()
    plt.plot(cdf_energies, cdf)
    plt.xscale('log')
    plt.xlabel('Energy [ev]')
    plt.title('Watt Spectrum CDF')
    plt.savefig('fission_spectrum_cdf.png')

    # Sample from fission CDF to get fission spectrum
    num_samples = 100000
    emitted_energies = np.zeros(num_samples)
    for i in range(num_samples):
        emitted_energies[i] = fissioner.emitNeutroneV()
    
    # Plot fission spectrum
    fig = plt.figure()
    plt.hist(emitted_energies, 100)
    plt.savefig('fission_spectrum.png')



    








