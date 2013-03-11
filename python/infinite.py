from pinspec import *
import numpy
import matplotlib.pyplot as matplt    # only need to import this for examples
import plotter
from SLBW import *

def main():


    # NOTE: If a user is going to homogenize materials
	#       then they must figure out the atom ratios 
	#       using geometric parameters and use that 
    #       when loading the isotope in a material

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 10000
    num_threads = 8
    log_setlevel(INFO)

    # Call SLBW to create XS
    nameoffile = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
    T=300 #Temp in Kelvin of target nucleus
    SLBWXS(nameoffile,T)

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')



    ############################################################################
    #EXAMPLE: How to plot thermal scattering CDFs 
    ############################################################################
    num_bins = h1.getNumThermalCDFBins()
    num_cdfs = h1.getNumThermalCDFs()
    Eprime_to_E = h1.retrieveEprimeToE(num_bins)
    E_to_kT = h1.retrieveEtokT(num_cdfs)

    cdfs = h1.retrieveThermalCDFs(num_cdfs*num_bins)
    cdfs = numpy.reshape(cdfs, [num_cdfs, num_bins])   # reshape 1D array to 2D
    dist = h1.retrieveThermalDistributions(num_cdfs*num_bins)
    dist = numpy.reshape(dist, [num_cdfs, num_bins])   # reshape 1D array to 2D

    # Plot the PDFs
    fig = matplt.figure()
    legend = []
    for i in range(num_cdfs):
        matplt.plot(Eprime_to_E, dist[i][:])
        legend.append(str(E_to_kT[i]) + ' kT')

    matplt.title(h1.getIsotopeName() + ' Thermal Scattering PDFs')
    matplt.ylabel('Probability')
    matplt.xlabel('Eprime / E')
    matplt.legend(legend)
    matplt.savefig(h1.getIsotopeName() + '_thermal_scattering_pdfs.png')


    # Plot the CDFs
    fig = matplt.figure()
    legend = []
    for i in range(num_cdfs):
        matplt.plot(Eprime_to_E, cdfs[i][:])
        legend.append(str(E_to_kT[i]) + ' kT')

    matplt.title(h1.getIsotopeName() + ' Thermal Scattering CDFs')
    matplt.ylabel('Cumulative Probability')
    matplt.xlabel('Eprime / E')
    matplt.legend(legend)
    matplt.savefig(h1.getIsotopeName() + 'thermal_scattering_cdfs.png')


    
    # Define materials
    mix = Material()
    mix.setMaterialName('fuel moderator mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(u238, 0.50)
    mix.addIsotope(u235, .025)
    
    log_printf(INFO, "Added isotopes")


    # Define regions
    region_mix = Region('infinite medium fuel-moderator mix', INFINITE)
    region_mix.setMaterial(mix)

    log_printf(INFO, "Made mixture region")
        
    # plot the fission spectrum the CDF
    plotter.plotFissionSpectrum()


    # Create a tally for the flux
    flux = Tally('total flux', REGION, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 2000, LOGARITHMIC)
    region_mix.addTally(flux)

    ############################################################################
    #EXAMPLE: How to set tally bin edges 
    ############################################################################
    # Create a tally for the absorption rate
    abs_rate = Tally('absorption rate', MATERIAL, ABSORPTION_RATE)
    abs_rate_bin_edges = numpy.array([0.1, 1., 5., 10., 100., 1000.])
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

    
    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics('Infinite_MC_Statistics', 'test')

    # Plotting
    plotter.plotFlux(flux)
    plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission', 'absorption'])
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission', 'absorption'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    plotter.plotMacroXS(mix, ['capture', 'elastic', 'fission', \
                                                    'absorption', 'total'])


if __name__ == '__main__':
    
    main()
