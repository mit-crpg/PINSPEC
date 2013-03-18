import matplotlib.pyplot as plt
import numpy as np
from pinspec import *
from process import *

# Series of functions use to plot XS data (energy, xs) that is saved
# in text files output by pinspec. The plotter creates
# a directory to store the text files, reads the text files
# and outputs png plots, and deletes the directory
# used to store the text files.

# Function to plot the microscopic cross section for a
# given isotope and array of reactions
def plotMicroXS(isotope, rxns, dir = '.'):
    
    # set input and output file names
    filename = dir + '/' + isotope.getIsotopeName() + '-micro-xs.png'
    
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
    plt.xlabel('Energy [eV]')
    plt.ylabel('$\sigma$' +' [b]')
    plt.title(isotope.getIsotopeName() + ' Microscopic XS')
    plt.legend(rxns, loc='lower left')
    plt.grid()
    fig.savefig(filename)


# function to plot the macroscopic cross section for a
# given material and arrary of reactions
def plotMacroXS(material, rxns, dir = '.'):
    
    # set input and output file names
    filename = dir + '/' + material.getMaterialName() + '-macro-xs.png'
    
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
    plt.xlabel('Energy [eV]')
    plt.ylabel('$\Sigma$'+' ['+'cm$^{-1}$'+']')
    plt.title(material.getMaterialName() + ' Macroscopic XS')
    plt.legend(rxns, loc='lower left')
    plt.grid()
    fig.savefig(filename)

# Function that plots the flux spectrum
def plotFlux(flux, dir = '.'):

    fig = plt.figure()

    num_bins = flux.getNumBins()
    flux_bin_centers = flux.retrieveTallyCenters(num_bins)
    flux_mu = flux.retrieveTallyMu(num_bins)

    # Plot the flux
    plt.plot(flux_bin_centers, flux_mu, lw=1)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux')
    plt.title('Batch-Averaged Flux')
    plt.grid()
    filename = dir + '/' + flux.getTallyName() + '_flux.png'
    fig.savefig(filename)


# Function that plots the flux spectrum for several tallies
def plotFluxes(fluxes, dir = '.'):

    fig = plt.figure()
    legend = []
    filename = dir + '/'

    for flux in fluxes:
        num_bins = flux.getNumBins()
        flux_bin_centers = flux.retrieveTallyCenters(num_bins)
        flux_mu = flux.retrieveTallyMu(num_bins)

        # Plot the flux
        plt.plot(flux_bin_centers, flux_mu, lw=1)
        legend.append(flux.getTallyName())
        filename += flux.getTallyName() + '_'

    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux')
    plt.title('Batch-Averaged Flux')
    plt.legend(legend, loc='lower right', prop={'size':12})
    plt.grid()
    fig.savefig(filename + 'flux.png')


# Function to plot the thermal scattering PDFs and CDFs
def plotThermalScattering(isotope, dir = '.'):
   
    num_bins = isotope.getNumThermalCDFBins()
    num_cdfs = isotope.getNumThermalCDFs()
    Eprime_to_E = isotope.retrieveEprimeToE(num_bins)
    E_to_kT = isotope.retrieveEtokT(num_cdfs)

    cdfs = isotope.retrieveThermalCDFs(num_cdfs*num_bins)
    cdfs = np.reshape(cdfs, [num_cdfs, num_bins])   # reshape 1D array to 2D
    dist = isotope.retrieveThermalDistributions(num_cdfs*num_bins)
    dist = np.reshape(dist, [num_cdfs, num_bins])   # reshape 1D array to 2D

    # Plot the PDFs
    fig = plt.figure()
    legend = []
    for i in range(num_cdfs):
        plt.plot(Eprime_to_E, dist[i][:])
        legend.append("%.2e" % E_to_kT[i] + ' kT')

    plt.title(isotope.getIsotopeName() + ' Thermal Scattering PDFs')
    plt.ylabel('Probability')
    plt.xlabel('E'+"'"+'/ E')
    plt.legend(legend, ncol=2, loc='upper right', prop={'size':14})
    plt.grid()
    filename = dir + '/' + isotope.getIsotopeName() + '_thermal_scattering_pdfs.png'
    plt.savefig(filename)

    # Plot the CDFs
    fig = plt.figure()
    legend = []
    for i in range(num_cdfs):
        plt.plot(Eprime_to_E, cdfs[i][:])
        legend.append("%.2e" % E_to_kT[i] + ' kT')
    plt.title(isotope.getIsotopeName() + ' Thermal Scattering CDFs')
    plt.ylabel('Cumulative Probability')
    plt.xlabel('E'+"'"+' / E')
    plt.legend(legend, ncol=2, loc='lower right', prop={'size':12})
    plt.grid()
    filename = dir + '/' + isotope.getIsotopeName() + 'thermal_scattering_cdfs.png'
    plt.savefig(filename)


# Function to plot the fission CDF and sample from the
# fission CDF to generate a fission spectrum
def plotFissionSpectrum(dir = '.'):
    
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
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Cumulative Probability')
    plt.title('Watt Spectrum CDF')
    plt.grid()
    filename = dir + '/' + 'fission_spectrum_cdf.png'
    plt.savefig(filename)

    # Sample from fission CDF to get fission spectrum
    num_samples = 10000000
    emitted_energies = np.zeros(num_samples)
    for i in range(num_samples):
        emitted_energies[i] = fissioner.emitNeutronMeV()

    # Bin the samples    
    binned_samples, bin_edges = np.histogram(emitted_energies, bins=1000, density=True)
    bin_centers = np.zeros(bin_edges.size-1)
    for i in range(bin_edges.size-1):
        bin_centers[i] = (bin_edges[i] + bin_edges[i+1]) / 2.0

    # Plot fission spectrum
    fig = plt.figure()
    plt.plot(bin_centers, binned_samples)
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Probability')
    plt.title('Watt Spectrum PDF')
    plt.grid()
    filename = dir + '/fission_spectrum_pdf.png'
    plt.savefig(filename)
    
    
def plotRI(RI, dir='.'):
    
    # Plot Resonance Integrals
    fig = plt.figure()
    bins = RI.bin_edges
    plt.semilogx(bins[0:-1], RI.RIs, drawstyle='steps')
    plt.xlabel('Energy [eV]')
    plt.ylabel('RI')
    plt.title('Resonance Integrals')
    plt.grid()
    filename = dir + '/RI.png'
    plt.savefig(filename)

def plotGroupXS(group_xs, dir='.'):
    
    # Plot Resonance Integrals
    fig = plt.figure()
    bins = group_xs.bin_edges
    plt.semilogx(bins[0:-1], group_xs.groupXS[:,:], drawstyle='steps')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Group XS')
    plt.title('Group XS')
    plt.grid()
    filename = dir + '/group_xs.png'
    plt.savefig(filename)    
    
    
    
    
    
    
    
    
    
    
