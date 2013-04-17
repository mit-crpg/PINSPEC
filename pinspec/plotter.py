##
# @file plotter.py
# @package pinspec.plotter
# @brief The plotter module provides utility functions to plot data from
#        PINSPEC's C++ classes, in particular, tally data and cross-sections,
#        thermal scattering PDFs/CDFs, etc.
#
# @author Samuel Shaner (shaner@mit.edu)
# @author William Boyd (wboyd@mit.edu
# @date March 10, 2013


import matplotlib.pyplot as plt
import numpy as np
from pinspec import *
from process import *
import os

## A static variable to auto-generate unique filenames for flux plots
flux_plot_num = 0

## A static variable for the output directory in which to save plots
subdirectory = "/plots/"
    

##
# @brief Plots the microscopoic cross-section(s) for one or more reaction
#        rates for an isotope.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotMicroXS(my_isotope, 'capture')
# @endcode
#
# @param isotope the isotope of interest
# @param rxns an array of reaction rate types (ie, ['capture, 'elastic'])
# @param loglog an optional argument to use a log-log scale
# @param uselegend an optional argument boolean to include a legend
# @param title an optional argument string with the plot title
# @param filename an optional argument string with the plot filename
def plotMicroXS(isotope, rxns, loglog=True, uselegend=True, \
                                            title='', filename=''):
        
    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)


    # make figure
    fig = plt.figure()
        
    # loop over rxns and plot
    num_energies = 0

    # Loop over all reaction rate types
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
        filename = directory + isotope.getIsotopeName() + '-micro-xs.png'
    else:
        filename = directory + filename.replace(' ', '-').lower() + '.png'

    fig.savefig(filename)



##
# @brief Plots the macroscopic cross-section(s) for one or more reaction
#        rates for a material.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotMacroXS(my_material, 'capture')
# @endcode
#
# @param material the material of interest
# @param rxns an array of reaction rate types (ie, ['capture, 'elastic'])
# @param loglog an optional argument to use a log-log scale
# @param uselegend an optional argument boolean to include a legend
# @param title an optional argument string with the plot title
# @param filename an optional argument string with the plot filename
def plotMacroXS(material, rxns, loglog=True, \
                                    uselegend=True, title='', filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory
        
    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

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
        filename = directory + material.getMaterialName() + '-macro-xs.png'
    else:
        filename = directory + filename.replace(' ', '-') + '.png'

    fig.savefig(filename)


##
# @brief Plots one or more flux tallies by energy.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotFlux(my_flux, title='one flux')
#         pinspec.plotter.plotFlux([flux1, flux2], title='two fluxes')
# @endcode
#
# @param fluxes the flux tall(ies) of interest
# @param loglog an optional argument boolean to use a log-log scale
# @param uselegend an optional argument boolean to include a legend
# @param title an optional argument string with the plot title
# @param filename an optional argument string with the plot filename
def plotFlux(fluxes, loglog=True, uselegend=False, title='', filename=''):

    global flux_plot_num
    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Make fluxes a list if it is not already
    if type(fluxes) is not list:
        fluxes = [fluxes]
    
    legend = []
    fig = plt.figure()


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
        filename = directory + '/' + filename.replace(' ', '-') + '.png'

    fig.savefig(filename)


##
# @brief Plots the thermal scattering PDFs and CDFs for an isotope.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotThermalScattering(my_isotope)
# @endcode
#
# @param isotope the isotope of interest
# @param uselegend an optional argument boolean to include a legend
# @param title an optional argument string with the plot title
# @param filename an optional argument string with the plot filename
def plotThermalScattering(isotope, uselegend=True, title='', filename=''):
   
    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

    num_bins = isotope.getNumThermalCDFBins()
    num_cdfs = isotope.getNumThermalCDFs()
    Eprime_to_E = isotope.retrieveEprimeToE(num_bins)
    E_to_kT = isotope.retrieveEtokT(num_cdfs)

    cdfs = isotope.retrieveThermalCDFs(num_cdfs*num_bins)
    cdfs = np.reshape(cdfs, [num_cdfs, num_bins])   # reshape 1D array to 2D
    dist = isotope.retrieveThermalPDFs(num_cdfs*num_bins)
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
        pdfsfilename = directory + isotope.getIsotopeName() + \
                                    '-thermal-scattering-pdfs.png'
    else:
        pdfsfilename = directory + filename.replace(' ', '-').lower() + \
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
        cdfsfilename = directory + isotope.getIsotopeName() + \
                                        '-thermal-scattering-cdfs.png'
    else:
        cdfsfilename = directory + filename.replace(' ', '-').lower() + \
                                        '-thermal-scattering-cdfs.png'

    plt.savefig(cdfsfilename)


##
# @brief Plots the fission spectrum PDF and CDF.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotFissionSpectrum()
# @endcode
def plotFissionSpectrum():
    
    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

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
    binned_samples, bin_edges = np.histogram(emitted_energies, bins=1000, 
                                             density=True)
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
    
    
##
# @brief Plots a resonance integral (RIEff or RITrue) as a step function.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotRI(my_resonance_integral)
# @endcode
def plotRI(RI, title='', filename=''):
    
    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

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
        filename = directory + 'RI.png'
    else:
        filename = directory + filename.replace(' ', '-').lower() +'-RI.png'

    plt.savefig(filename)


##
# @brief Plots a multi-group cross-section as a step function.
# @details This method generates and saves the plot as a *.png in the
#          plotting output directory. A user may invoke this function from
#          a PINSPEC Python file as follows:
#
# @code
#         pinspec.plotter.plotGroupXS(my_group_xs)
# @endcode
def plotGroupXS(group_xs, title='', filename=''):
    
    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

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
        filename = directory + filename.replace(' ', '-').lower() + \
                                                        '-group-xs.png'

    plt.savefig(filename)



##
# @brief This method takes in a region or a geometry object and plots a
#        color-coded 2D surface plot representing a planar slice through 
#         the space.
# @param space a bounded region or heterogeneous geometry object
# @param plane the optional 'xy', 'xz' or 'yz' plane the user wishes to plot
# @param loc the optional x, y, or z position of the planar plot
# @param lim1 the optional bounding limits for the first planar dimension
# @param lim2 the optional bounding limits for the second planar dimension
# @param gridsize an optional number of grid cells for the plot
# @param filename an optional filename
def plotSlice(space, plane='XY', loc=0.0, lim1=[-2., 2.], lim2=[-2., 2.], 
              gridsize=100, filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

    # Error checking
    if not isinstance(space, (BoundedModeratorRegion, BoundedFuelRegion, \
                                   BoundedGeneralRegion, Geometry)):
        py_printf('ERROR', 'Unable to plot a %s slice since input ' + \
                  'did not contain a region or geometry', plane)
    if not isinstance(plane, str):
        py_printf('ERROR', 'Unable to plot a %s slice for space %s.' + \
                      'PINSPEC only supports plots of XY, XZ and YZ slices.',
                      str(plane), space.getName())
    if plane.lower() != 'xy' and plane.lower() != 'xz' \
            and plane.lower() != 'yz':
        print 'plane.lower = ' + str(plane.lower())
        py_printf('ERROR', 'Unable to plot a %s slice for space %s.' + \
                      'PINSPEC only supports plots of XY, XZ and YZ slices.',
                      plane, space.getName()) 
    if len(lim1) is not 2:
        py_printf('ERROR', 'Unable to plot a %s slice for space %s since ' + \
                      'the lim1 parameter is of length %d rather than 2', 
                      plane, space.getName(), len(lim1))
    if len(lim2) is not 2:
        py_printf('ERROR', 'Unable to plot a %s slice for space %s since ' + \
                      'the lim2 parameter is of length %d rather than 2', 
                      plane, space.getName(), len(lim2))
    if not isinstance(lim1[0], float) or not isinstance(lim1[1], float):
        py_printf('ERROR', 'Unable to plot a %s slice for space %s since ' + \
                      'lim1 does not contain floating point values', plane,
                      space.getName())
    if not isinstance(lim2[0], float) or not isinstance(lim2[1], float):
        py_printf('ERROR', 'Unable to plot a %s slice for space %s since ' + \
                      'lim2 does not contain floating point values', plane,
                      space.getName())
    if lim1[0] >= lim1[1]:
        py_printf('ERROR', 'Unable to plot a %s slice for space %s since ' + \
                      'the minimum lim1 value (%f) is greater than the ' + \
                      'maximum lim1 value (%f)', space.getName(), plane,
                      lim1[0], lim1[1])
    if lim2[0] >= lim2[1]:
        py_printf('ERROR', 'Unable to plot a %s slice for space %s since ' + \
                      'the minimum lim2 value (%f) is greater than the ' + \
                      'maximum lim2 value (%f)', space.getName(), plane,
                      lim2[0], lim2[1])
    if not isinstance(gridsize, int):
        py_printf('ERROR', 'Unable to plot a %s slice for space %s ' + \
                      'Since the gridsize %s is not an integer', plane,
                      space.getName(), str(gridsize))
    if gridsize <= 0:
        py_printf('Error', 'Unable to plot a %s slice for space %s ' + \
                      'with a negative gridsize (%d)', plane,
                      space.getName(), gridsize)
    if not isinstance(loc, float):
        py_printf('Error', 'Unable to plot a %s slice for space %s ' + \
                      'since the input location of the slice is not a number',
                  plane, space.getName())

    # Initialize a numpy array for the surface colors
    surface = numpy.zeros((gridsize, gridsize), dtype=np.int32)

    dim1 = np.linspace(lim1[0], lim1[1], gridsize)
    dim2 = np.linspace(lim2[0], lim2[1], gridsize)

    if plane.lower() == 'xy':

        # We are plotting a region
        if isinstance(space, (BoundedRegion, BoundedModeratorRegion, \
                                  BoundedFuelRegion)):
            for i in range(len(dim1)):
                for j in range(len(dim2)):
            
                    # Color space for points within the space
                    if space.contains(dim1[i], dim2[j], loc):
                        surface[j][i] = 1.0

        # We are plotting the geometry
        else:

            for i in range(len(dim1)):
                for j in range(len(dim2)):
            
                    # Color space for points within the space
                    if space.contains(dim1[i], dim2[j], loc):
                        region = space.findContainingRegion(dim1[i], \
                                                            dim2[j], loc)
                        surface[j][i] = region.getUid()

    if plane.lower() == 'xz':

        # We are plotting a region
        if isinstance(space, (BoundedRegion, BoundedModeratorRegion, \
                                  BoundedFuelRegion)):
            for i in range(len(dim1)):
                for j in range(len(dim2)):
            
                    # Color space for points within the space
                    if space.contains(dim1[i], loc, dim2[j]):
                        surface[j][i] = 1.0

        # We are plotting the geometry
        else:

            for i in range(len(dim1)):
                for j in range(len(dim2)):
            
                    # Color space for points within the space
                    if space.contains(dim1[i], loc, dim2[j]):
                        region = space.findContainingRegion(dim1[i], loc,
                                                            dim2[j])
                        surface[j][i] = region.getUid()


    if plane.lower() == 'yz':

        # We are plotting a region
        if isinstance(space, (BoundedRegion, BoundedModeratorRegion, \
                                  BoundedFuelRegion)):

            # We are plotting a region
            for i in range(len(dim1)):
                for j in range(len(dim2)):
            
                    # Color space for points within the space
                    if space.contains(loc, dim1[i], dim2[j]):
                        surface[j][i] = 1.0

        # We are plotting the geometry
        else:

            for i in range(len(dim1)):
                for j in range(len(dim2)):
            
                    # Color space for points within the space
                    if space.contains(loc, dim1[i], dim2[j]):
                        region = space.findContainingRegion(loc, dim1[i],
                                                            dim2[j])
                        surface[j][i] = region.getUid()

    
    fig = plt.figure()
    plt.pcolor(dim1, dim2, surface)
    plt.axis([lim1[0], lim1[1], lim2[0], lim2[1]])

    if plane.lower() == 'xy':
        plt.title('' + plane.upper() + ' Plane at z = ' + str(loc))
    if plane.lower() == 'xz':
        plt.title('' + plane.upper() + ' Plane at y = ' + str(loc))
    else:
        plt.title('' + plane.upper() + ' Plane at x = ' + str(loc))

    if filename is '':
        filename = directory + '/' + plane.lower() + '-slice.png'
    else:
        filename = directory + filename.replace(' ', '-').lower() + \
                      '-' + plane.lower() + '-slice.png'

    plt.savefig(filename)


##
# @brief Plots a vector field of the fission site source distribution
# @param geometry a HETEROGENEOUS geometry class object
# @param num_samples an optional number of samples
# @param filename an optional filename
def plotFissionSourceDist(geometry, num_samples=1000, filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

    # Error checking
    if not isinstance(geometry, Geometry):
        py_printf('ERROR', 'Unable to plot the fission source distribution ' + \
                  'since the parameters did not include a geometry class object')
    if geometry.getSpatialType() != HETEROGENEOUS:
        py_printf('ERROR', 'Unable to plot the fission source distribution ' + \
                      'since for a non HETEROGENEOUS type geometry')
    if not isinstance(num_samples, int):
        py_printf('ERROR', 'Unable to plot the fission source distribution ' + \
                      'since a non-integer number of samples was input')
    if num_samples < 0:
        py_printf('ERROR', 'Unable to plot the fission source distribution ' + \
                   'for %d number of particle since it is negative', num_samples)

    neutron = createNewNeutron()

    x = []
    y = []
    u = []
    v = []


    # Sample fission source neutrons and save their locations and velocity vectors
    for i in range(1000):
        geometry.initializeSourceNeutron(neutron)
    
        x.append(neutron._x)
        y.append(neutron._y)
        u.append(neutron._u)
        v.append(neutron._v)


    # Plot the vector field
    plt.figure()
    plt.quiver(x,y,u,v,angles='xy', color='r')
    plt.title('Fission Source Distribution')

    if filename is '':
        filename = directory + '/' + 'fission-site-dist.png'
    else:
        filename = directory + filename.replace(' ', '-').lower() + '-.png'

    plt.savefig(filename)


##
# @brief This method takes in a geometry class object and tracks a neutron
#        as it travels across the geomgetry.
# @details A color-coded 2D surface plot representing a planar slice through 
#         the space is produced as the backdrop to a series of line segments
#         representing a neutron's path across the pin cell.
# @param geometry a heterogeneous geometry object
# @param plane the optional 'xy', 'xz' or 'yz' plane the user wishes to plot
# @param num_moves the number of neutron moves across the geometry to plot
# @param loc the optional x, y, or z position of the planar plot
# @param lim1 the optional bounding limits for the first planar dimension
# @param lim2 the optional bounding limits for the second planar dimension
# @param gridsize an optional number of grid cells for the plot
# @param filename an optional filename
def trackANeutron(geometry, plane='XY', num_moves=100, loc=0.0, \
                      lim1=[-2., 2.], lim2=[-2., 2.], gridsize=100, filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
            os.makedirs(directory)

    # Error checking
    if not isinstance(geometry, Geometry):
        py_printf('ERROR', 'Unable to track a neutron since input ' + \
                      'did not contain a geometry class object')
    if geometry.getSpatialType() != HETEROGENEOUS:
        py_printf('ERROR', 'Unable to track a neutron since geometry is ' + \
                      'not a HETEROGENEOUS type geometry')
    if not isinstance(num_moves, int):
        py_printf('ERROR', 'Unable to track a neutrons since the number ' + \
                      'of moves input is %s which is not an integer', \
                      str(num_moves))
    if (num_moves < 0):
        py_printf('ERROR', 'Unable to track a neutron since the number ' + \
                      'of moves input is %d which is negative', num_moves)
    if not isinstance(plane, str):
        py_printf('ERROR', 'Unable to track a neutron since for slice %s ' + \
                      'PINSPEC only supports plots of XY, XZ and YZ slices.',
                      str(plane))
    if plane.lower() != 'xy' and plane.lower() != 'xz' \
            and plane.lower() != 'yz':
        print 'plane.lower = ' + str(plane.lower())
        py_printf('ERROR', 'Unable to track a neutron since for slice %s ' + \
                      'PINSPEC only supports plots of XY, XZ and YZ slices.',
                      plane) 
    if len(lim1) is not 2:
        py_printf('ERROR', 'Unable to track a neutron since for slice %s ' + \
                      'since the lim1 parameter is of length %d rather than 2', 
                      plane, len(lim1))
    if len(lim2) is not 2:
        py_printf('ERROR', 'Unable to  track a neutron since for slice %s ' + \
                      'since the lim2 parameter is of length %d rather than 2', 
                      plane,  len(lim2))
    if not isinstance(lim1[0], float) or not isinstance(lim1[1], float):
        py_printf('ERROR', 'Unable to  track a neutron since for slice %s ' + \
                      'since lim1 does not contain floating point values', plane)
    if not isinstance(lim2[0], float) or not isinstance(lim2[1], float):
        py_printf('ERROR', 'Unable to track a neutron since for slice %s ' + \
                      'since lim2 does not contain floating point values', plane)
    if lim1[0] >= lim1[1]:
        py_printf('ERROR', 'Unable to track a neutron since for slice %s ' + \
                      'since the minimum lim1 value (%f) is greater than the ' + \
                      'maximum lim1 value (%f)', plane, lim1[0], lim1[1])
    if lim2[0] >= lim2[1]:
        py_printf('ERROR', 'Unable to track a neutron since for slice %s ' + \
                      'since the minimum lim2 value (%f) is greater than the ' + \
                      'maximum lim2 value (%f)', plane, lim2[0], lim2[1])
    if not isinstance(gridsize, int):
        py_printf('ERROR', 'Unable to  track a neutron since for slice %s ' + \
                      'since the gridsize %s is not an integer', plane,
                      str(gridsize))
    if gridsize <= 0:
        py_printf('Error', 'Unable to  track a neutron since for slice %s ' + \
                      'for a negative gridsize (%d)', plane, gridsize)
    if not isinstance(loc, float):
        py_printf('Error', 'Unable to track a neutron for slice %s ' + \
                  'since the input location of the slice is not a number', plane)


    # First plot the regions in the geometry as a backdrop for the
    # neutron's path

    # Initialize a numpy array for the surface colors
    surface = numpy.zeros((gridsize, gridsize), dtype=np.int32)

    # Initialize 
    dim1 = np.linspace(lim1[0], lim1[1], gridsize)
    dim2 = np.linspace(lim2[0], lim2[1], gridsize)

    if plane.lower() == 'xy':

        for i in range(len(dim1)):
            for j in range(len(dim2)):
            
                # Color space for points within the space
                if geometry.contains(dim1[i], dim2[j], loc):
                    region = geometry.findContainingRegion(dim1[i], \
                                                        dim2[j], loc)
                    surface[j][i] = region.getUid()

    if plane.lower() == 'xz':

        for i in range(len(dim1)):
            for j in range(len(dim2)):
            
                # Color space for points within the space
                if geometry.contains(dim1[i], loc, dim2[j]):
                    region = geometry.findContainingRegion(dim1[i], loc,
                                                        dim2[j])
                    surface[j][i] = region.getUid()


    if plane.lower() == 'yz':

        for i in range(len(dim1)):
            for j in range(len(dim2)):
         
                # Color space for points within the space
                if geometry.contains(loc, dim1[i], dim2[j]):
                    region = geometry.findContainingRegion(loc, dim1[i],
                                                        dim2[j])
                    surface[j][i] = region.getUid()

    
    # Plot the regions in the geometry
    fig = plt.figure()
    plt.pcolor(dim1, dim2, surface)
    plt.axis([lim1[0], lim1[1], lim2[0], lim2[1]])

    # Create a neutron to track
    neutron = createNewNeutron()
    geometry.initializeSourceNeutron(neutron)

    # Initialize empty lists to track the neutron's location
    x = []
    y = []
    z = []

    # Track the neutron
    for i in range(num_moves):

        # If the neutron was killed on the last collision, then break loop
        if not neutron._alive:
            break

        geometry.findContainingRegion(neutron)
        region = neutron._region
        region.collideNeutron(neutron)
        x.append(neutron._x)
        y.append(neutron._y)
        z.append(neutron._z)

    # Plot the neutron path throughout the geometry
    if plane.lower() == 'xy':
        plt.plot(x,y,'b-')
        plt.plot(x,y,'bo', markersize=7, markeredgecolor='w')
    elif plane.lower() == 'yz':
        plt.plot(y,z,'b-') 
        plt.plot(y,z,'bo', markersize=7, markeredgecolor='w')
    else:
        plt.plot(x,z,'b-') 
        plt.plot(x,z,'bo', markersize=7, markeredgecolor='w')

    plt.title('Tracking a Neutron in ' + plane.upper())

    if filename is '':
        filename = directory + '/' + 'track-a-neutron-' + plane.lower() + '.png'
    else:
        filename = directory + filename.replace(' ', '-').lower() + \
                         '-' + plane.lower() + '-.png'

    plt.savefig(filename)
