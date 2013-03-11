import matplotlib.pyplot as plt
from numpy import *
from pinspec import *

# Series of functions use to plot XS data (energy, xs) that is saved
# in text files output by pinspec. The plotter creates
# a directory to store the text files, reads the text files
# and outputs png plots, and deletes the directory
# used to store the text files.

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











