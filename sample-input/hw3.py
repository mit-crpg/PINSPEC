import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

def main():

    ###########################################################################
    ###################   Initialize Params and Isotopes   ####################
    ###########################################################################

    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 100000
    num_threads = 4
    setOutputDirectory('HW3');
    log_setlevel(INFO)

    py_printf('TITLE', 'Simulation of homework 3 for 2012 22.211')
    py_printf('INFO', 'Initializing isotopes...')

	# Initialize isotopes
    h1 = Isotope('H-1')
    u238 = Isotope('U-238')

    # Create an artifical capture xs for hydrogen
    norm_const = 0.025
    h1_capture_energies = numpy.logspace(-5., 7.5, 500);
    h1_capture_xs = numpy.sqrt(norm_const/h1_capture_energies) * 2.;
    h1.setCaptureXS(h1_capture_energies, h1_capture_xs)

	# Zero out the scatter, fission xs for U-238
    energy = numpy.array([1E-7, 2e7])	# energy bounds
    xs = numpy.array([0.0])				# one group xs
    u238.setMultigroupElasticXS(energy, xs)
    u238.setMultigroupFissionXS(energy, xs)

    py_printf('INFO', 'Plotting microscopic cross-sections...')

    plotter.plotMicroXS(h1, ['capture', 'elastic', 'fission'])
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'])

	# Turn off thermal scattering for U-238
    u238.neglectThermalScattering()

    py_printf('INFO', 'Creating a fuel-moderator mixture material...')

    # Define materials
    mix = Material('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(u238, 1E-1)

    py_printf('INFO', 'Creating an infinite homogeneous region...')
    
    # Define regions
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    py_printf('INFO', 'Creating the geometry...')
 
    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    py_printf('INFO', 'Initializing tallies for flux and absorption rates...')

    # Create a tally for the flux
    flux1 = TallyFactory.createTally(geometry, FLUX, 'U/H = 1E-6')
    flux2 = TallyFactory.createTally(geometry, FLUX, 'U/H = 0.01')
    flux3 = TallyFactory.createTally(geometry, FLUX, 'U/H = 0.1')
    flux4 = TallyFactory.createTally(geometry, FLUX, 'U/H = 1.0')
    flux1.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
    flux2.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
    flux3.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
    flux4.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)
    fluxes = [flux1, flux2, flux3, flux4]

    # Create tallies to compute absorption in u238 and the material
    u238_abs_rate = TallyFactory.createTally(u238, ABSORPTION_RATE, 'u238 abs')
    tot_abs_rate = TallyFactory.createTally(region_mix, ABSORPTION_RATE, 'tot abs')
    abs_rate_flux = TallyFactory.createTally(region_mix, FLUX)

    abs_rate_bin_edges = numpy.array([1E-5, 1., 6., 10., 25., 50., 100., 1000.])
    tot_abs_rate.generateBinEdges(1E-7, 2E7, 1, EQUAL)
    u238_abs_rate.setBinEdges(abs_rate_bin_edges)
    abs_rate_flux.setBinEdges(abs_rate_bin_edges)

    py_printf('INFO', 'Registering tallies with the TallyBank...')

    TallyBank.registerTally(u238_abs_rate, region_mix)
    TallyBank.registerTally(tot_abs_rate, region_mix)
    TallyBank.registerTally(abs_rate_flux)


    ###########################################################################
    #########################   Problems 1 and 2   ############################
    ###########################################################################

    py_printf('HEADER', 'Problems 1 and 2')

    u_h_ratios = [1E-6, 1E-2, 1E-1, 1E0]
    temps = [300, 600, 900, 1200]
    abs_rate_ratios = []
    Eff_RIs = []
    True_RIs = []


    for temp in range(4):

        py_printf('TITLE', 'Temperature = %dK', temps[temp])

        Eff_RIs.append([])
        abs_rate_ratios.append([])

        SLBW.SLBWXS('U-238', temps[temp], 'capture')
        u238.loadXS('capture')
        RI = process.RITrue(u238, abs_rate_bin_edges, reaction='capture')
        RI.setName('True RI (Temp=%dK)' % temps[temp])
        True_RIs.append(RI)

        for u_h_ratio in range(4):

            TallyBank.registerTally(fluxes[u_h_ratio])

            # Update U-238 atom ratio
            mix.addIsotope(u238, u_h_ratios[u_h_ratio])

            # Run Monte Carlo simulation
            geometry.runMonteCarloSimulation();

            # Compute the resonance integrals
            RI = process.RIEff(abs_rate_flux, u238_abs_rate)
            RI.setName('Effective RI (U/H=%1.0E)' % u_h_ratios[u_h_ratio])
            Eff_RIs[temp].append(RI)

	        # Compute absorption rate ratios
            abs_rate_ratio = u238_abs_rate / tot_abs_rate
            abs_rate_ratio.setTallyName('U238 / Tot for U/H=%1.0E' \
                                                    % u_h_ratios[u_h_ratio])
            abs_rate_ratios[temp].append(abs_rate_ratio)

            # Append flux for this simulation to an array and deregister it
            # for next simulation
            fluxes[u_h_ratio].normalizeBatchMu()
            TallyBank.deregisterTally(fluxes[u_h_ratio])

        # Plot fluxes
        plotter.plotFluxes(fluxes,title='Flux for Temp = ' + str(temps[temp]) \
                            + 'K', filename='flux-temp-'+str(temps[temp])+'K')

        # print the reaction rate ratios and resonance integrals to the shell
        process.printTallies(abs_rate_ratios[temp], 
                                                header='U238 Abs. / Tot. Abs.')
        process.printRIs(Eff_RIs[temp], header='Temp=%dK'%temps[temp])

    process.printRIs(True_RIs, header='T=300-1200K')


    py_printf('TITLE', 'Finished')
    

if __name__ == '__main__':

    main()  

