import matplotlib.pyplot as plt
import numpy as np
from pinspec import *
from process import *
import os

# Series of functions use to plot XS data (energy, xs), thermal
# scattering PDFs and CDFs, flux, etc

flux_plot_num = 0

def createDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    

# Function to plot the microscopic cross section for a
# given isotope and array of reactions
def plotMicroXS(isotope, rxns, directory = '.', loglog=True, \
                    uselegend=True, title='', filename=''):

    createDirectory(directory)
        
    # make figure
    fig = plt.figure()
        
    # loop over rxns and plot
    num_energies = 0

    for rxn in rxns:

	    # retrieve xs and xsenergies
        num_energies = isotope.getNumXSEnergies(rxn)
        energies = isotope.retrieveXSEnergies(num_energies, rxn)
        xs = isotope.retrieveXS(num_energies, rxn)

        # plot xs
        plt.plot(energies, xs, lw=1)

    plt.xlabel('Energy [eV]')
    plt.ylabel('$\sigma$' +' [b]')
    plt.grid()

    if (uselegend):
        plt.legend(rxns, loc='lower left')

    if (loglog):
        plt.xscale('log')
        plt.yscale('log')

    plt.title(isotope.getIsotopeName() + ' Microscopic XS')

    if title is '':
        plt.title(isotope.getIsotopeName() + ' Microscopic XS')
    else:
        plt.title(title.title())

    if filename is '':
        filename = directory + '/' + isotope.getIsotopeName() + '-micro-xs.png'
    else:
        filename = directory +'/' + filename.replace(' ', '-').lower() + '.png'

    fig.savefig(filename)



# function to plot the macroscopic cross section for a
# given material and arrary of reactions
def plotMacroXS(material, rxns, directory = '.', loglog=True, \
                        uselegend=True, title='', filename=''):
    
    createDirectory(directory)
    
    # make figure
    fig = plt.figure()
    
    # loop over rxns and plot
    num_energies = 0
    for rxn in rxns:

	    # retrieve xs and xs energies
        num_energies = material.getNumXSEnergies(rxn)
        energies = material.retrieveXSEnergies(num_energies, rxn)
        xs = material.retrieveXS(num_energies, rxn)

        # plot xs
        plt.plot(energies, xs, lw=1)


    plt.xlabel('Energy [eV]')
    plt.ylabel('$\Sigma$'+' ['+'cm$^{-1}$'+']')
    plt.grid()
 
    if (uselegend):
        plt.legend(rxns, loc='lower left')

    if (loglog):
        plt.xscale('log')
        plt.yscale('log')

    if title is '':
        plt.title(material.getMaterialName() + ' Macroscopic XS')
    else:
        plt.title(title.title())

    if filename is '':
        filename = directory+'/' + material.getMaterialName() + '-macro-xs.png'
    else:
        filename = directory + '/' + filename.replace(' ', '-') + '.png'

    fig.savefig(filename)



# Function that plots the flux spectrum
def plotFlux(flux, directory = '.', loglog=True, title='', filename=''):

    global flux_plot_num

    if not flux.hasComputedBatchStatistics():
        flux.computeBatchStatistics()

    createDirectory(directory)

    fig = plt.figure()
        
    num_bins = flux.getNumBins()
    flux_bin_centers = flux.retrieveTallyCenters(num_bins)
    flux_mu = flux.retrieveTallyMu(num_bins)

    # Plot the flux
    plt.plot(flux_bin_centers, flux_mu, lw=1)

    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux')
    plt.grid()

    if (loglog):
        plt.xscale('log')
        plt.yscale('log')

    if title is '':
        plt.title('Batch-Averaged Flux')
    else:
        plt.title(title.title())

    if filename is '':
        filename = directory + '/' + 'flux-' + str(flux_plot_num) + '.png'
        flux_plot_num += 1
    else:
        filename = directory + '/' + filename.replace(' ', '-') + '.png'

    fig.savefig(filename)


# Function that plots the flux spectrum for several tallies
def plotFluxes(fluxes, directory = '.', loglog=True, uselegend=True, \
                                                    filename='', title=''):

    global flux_plot_num

    createDirectory(directory)

    fig = plt.figure()
    legend = []

    for flux in fluxes:

        if not flux.hasComputedBatchStatistics():
            flux.computeBatchStatistics()

        num_bins = flux.getNumBins()
        flux_bin_centers = flux.retrieveTallyCenters(num_bins)
        flux_mu = flux.retrieveTallyMu(num_bins)

        # Plot the flux
        plt.plot(flux_bin_centers, flux_mu, lw=1)
        legend.append(flux.getTallyName())


    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux')
    plt.grid()

    if (loglog):
        plt.xscale('log')
        plt.yscale('log')

    if (uselegend):
        plt.legend(legend, loc='lower right', prop={'size':12})

    if title is '':
        plt.title('Batch-Averaged Flux')
    else:
        plt.title(title.title())

    if filename is '':
        filename = directory + '/' + 'flux-' + str(flux_plot_num) + '.png'
        flux_plot_num += 1
    else:
        filename = directory + '/'+ filename.replace(' ', '-').lower() + '.png'

    fig.savefig(filename)


# Function to plot the thermal scattering PDFs and CDFs
def plotThermalScattering(isotope, directory = '.', uselegend=True, \
                                                title='', filename=''):
   
    createDirectory(directory)

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

    plt.ylabel('Probability')
    plt.xlabel('E'+"'"+'/ E')
    plt.grid()

    if (uselegend):
        plt.legend(legend, ncol=2, loc='upper right', prop={'size':14})

    if title is '':
        plt.title(isotope.getIsotopeName() + ' Thermal Scattering PDFs')
    else:
        plt.title(title.title() + ' Thermal Scattering PDFs')

    if filename is '':
        pdfsfilename = directory + '/' + isotope.getIsotopeName() + \
                                    '-thermal-scattering-pdfs.png'
    else:
        pdfsfilename = directory + '/' + filename.replace(' ', '-').lower() + \
                                    '-thermal-scattering-pdfs.png'

    plt.savefig(pdfsfilename)

    # Plot the CDFs
    fig = plt.figure()
    legend = []
    for i in range(num_cdfs):
        plt.plot(Eprime_to_E, cdfs[i][:])
        legend.append("%.2e" % E_to_kT[i] + ' kT')

    plt.ylabel('Cumulative Probability')
    plt.xlabel('E'+"'"+' / E')
    plt.grid()

    if (uselegend):
        plt.legend(legend, ncol=2, loc='lower right', prop={'size':12})

    if title is '':
        plt.title(isotope.getIsotopeName() + ' Thermal Scattering CDFs')
    else:
        plt.title(title.title() + ' Thermal Scattering CDFs')

    if filename is '':
        cdfsfilename = directory + '/' + isotope.getIsotopeName() + \
                                        '-thermal-scattering-cdfs.png'
    else:
        cdfsfilename = directory + '/' + filename.replace(' ', '-').lower() + \
                                        '-thermal-scattering-cdfs.png'

    plt.savefig(cdfsfilename)


# Function to plot the fission CDF and sample from the
# fission CDF to generate a fission spectrum
def plotFissionSpectrum(directory = '.'):
    
    createDirectory(directory)

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
    plt.grid()
    plt.title('Watt Spectrum CDF')
    filename = directory + '/' + 'fission-spectrum-cdf.png'
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
    plt.grid()
    plt.title('Watt Spectrum PDF')
    filename = directory + '/fission_spectrum_pdf.png'
    plt.savefig(filename)
    
    
def plotRI(RI, directory='.', title='', filename=''):
    
    createDirectory(directory)

    # Plot Resonance Integrals
    fig = plt.figure()
    bins = RI.bin_edges
    plt.semilogx(bins[0:-1], RI.RIs, drawstyle='steps-post')
    plt.xlabel('Energy [eV]')
    plt.ylabel('RI')
    plt.grid()

    if title is '':
        plt.title('Resonance Integrals')
    else:
        plt.title(title.title() + ' Resonance Integrals')

    if filename is '':
        filename = directory + '/RI.png'
    else:
        filename = directory+'/'+filename.replace(' ', '-').lower() +'-RI.png'

    plt.savefig(filename)


def plotGroupXS(group_xs, directory='.', title='', filename=''):
    
    createDirectory(directory)

    # Plot Resonance Integrals
    fig = plt.figure()
    bins = group_xs.bin_edges
    plt.semilogx(bins[0:-1], group_xs.groupXS[:,:], drawstyle='steps-post')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Group XS')
    plt.grid()

    if title is '':
        plt.title('Group XS')
    else:
        plt.title(title.title() + ' Group XS')

    if filename is '':
        filename = directory + '/group-xs.png'
    else:
        filename = directory+'/'+filename.replace(' ', '-').lower() + \
                                                        '-group-xs.png'

    plt.savefig(filename)
