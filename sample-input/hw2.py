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
    num_neutrons = 1000000
    setOutputDirectory('HW2');
    log_setlevel(INFO)

    py_printf('TITLE', 'Simulation of homework 2 for 2012 22.211')
    py_printf('INFO', 'Initializing isotopes...')

	# Initialize isotopes
    h1 = Isotope('H-1')
    c12 = Isotope('C-12')

    # Create an artifical capture xs for hydrogen
    norm_const = 0.025
    h1_capture_energies = numpy.logspace(-5., 7.5, 500);
    h1_capture_xs = numpy.sqrt(norm_const/h1_capture_energies) * 7.;
    h1.setCaptureXS(h1_capture_energies, h1_capture_xs)

    py_printf('INFO', 'Plotting microscopic cross-sections...')
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'fission'])
    plotter.plotMicroXS(c12, ['capture', 'elastic', 'fission'])

    ###########################################################################
    ##########################   Problems 1 and 2   ###########################
    ###########################################################################

    py_printf('HEADER', 'Problems 1 and 2')
    py_printf('INFO', 'Plotting thermal scattering distributions...')

    #Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScattering(h1)

    #Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScattering(c12)

    ###########################################################################
    ###########################   Problems 3-6   ##############################
    ###########################################################################

    py_printf('HEADER', 'Problems 3-6')
    py_printf('INFO', 'Initializing the collision rate tally...')

    h1_material = Material('H-1')
    h1_material.setDensity(0.07778, 'g/cc')
    h1_material.addIsotope(h1, 2.0)
    
    py_printf('INFO', 'Initializing tallies for the flux, collision rate, etc')

    flux = TallyFactory.createTally(h1_material, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 1000, LOGARITHMIC)

    coll_rate_1eV = TallyFactory.createTally(h1_material, COLLISION_RATE)
    coll_rate_1eV.generateBinEdges(1E-1, 2E6, 1, EQUAL)

    coll_rate = TallyFactory.createTally(h1_material, COLLISION_RATE)
    coll_rate.generateBinEdges(1E-7, 2E6, 1, EQUAL)

    times = TallyFactory.createTally(h1_material, INTERCOLLISION_TIME)
    times.generateBinEdges(1E-7, 2E6, 1, EQUAL)

    py_printf('INFO', 'Simulating %d neutrons in H-1...', num_neutrons)

    fissioner = Fissioner()
    neutron = initializeNewNeutron()

    for i in range(num_neutrons):

        # Sample a fission energy from the Watt fission spectrum 
        neutron._energy = fissioner.emitNeutroneV()
        neutron._alive = True
        reached_one_ev = False

        # Simulate neutron until it is absorbed in H-1
        while(neutron._alive):

            h1_material.collideNeutron(neutron)
            flux.tally(neutron)
            times.tally(neutron)
            coll_rate.tally(neutron)

            if neutron._energy < 1.0:
                reached_one_ev = True

            if not reached_one_ev:
                coll_rate_1eV.tally(neutron)


    py_printf('INFO', 'Plotting the flux...')

    plotter.plotFlux(flux, title='H-1 Flux', filename='h-1-flux')

    num_collisions = process.computeMeanNumCollisions(coll_rate_1eV, \
                                                              num_neutrons)

    py_printf('RESULT', 'Mean # of collisions to 1 eV: %f', num_collisions)

    num_collisions = process.computeMeanNumCollisions(coll_rate, \
                                                              num_neutrons)

    py_printf('RESULT', 'Mean # of collisions to death: %f', num_collisions)

    mean_lifetime = process.computeMeanNeutronLifetime(times, num_neutrons)

    py_printf('RESULT', 'Avg neutron lifetime: %1.2E seconds', mean_lifetime)

    py_printf('TITLE', 'Finished')


if __name__ == '__main__':

    main()  

