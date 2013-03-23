import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

def main():

    ###########################################################################
    ###################   Initialize Params and Isotopes   ###################
    ###########################################################################

    # Set main simulation params
    num_neutrons = 10000
    output_dir = 'HW1'
    log_setlevel(INFO)

    py_printf('INFO', 'Initializing isotopes...')

	# Initialize isotopes
    h1 = Isotope('H-1')
    c12 = Isotope('C-12')

	# Zero out the capture xs
    energy = numpy.array([1E-7, 2e7])	# energy bounds
    xs = numpy.array([0.0])				# one group xs
    h1.setMultigroupCaptureXS(energy, xs)
    c12.setMultigroupCaptureXS(energy, xs)

    py_printf('INFO', 'Plotting microscopic cross-sections...')

    plotter.plotMicroXS(h1, ['capture', 'elastic', 'fission'], output_dir )
    plotter.plotMicroXS(c12, ['capture', 'elastic', 'fission'], output_dir)

	# Turn off thermal scattering
    h1.neglectThermalScattering()
    c12.neglectThermalScattering()


    ###########################################################################
    #########################   Problems 1 and 2   ############################
    ###########################################################################

    num_generations = 15
    max_energy = 2E6

    py_printf('NORMAL', 'Beginning Problems 1 and 2...')
    py_printf('NORMAL', '# generations = %d\t\t# neutrons / generation = %d', \
                                                num_generations, num_neutrons)
    py_printf('INFO', 'Initializing generational flux tallies...')

	# Create flux tallies for each generation
    h1_fluxes = []
    c12_fluxes = []

    for i in range(num_generations):

        h1_fluxes.append(TallyFactory.createTally('h1 flux', h1, COLLISION_RATE))
        h1_fluxes[i].generateBinEdges(0., 2E6, 1000, EQUAL)

        c12_fluxes.append(TallyFactory.createTally('c12 flux', c12, COLLISION_RATE))
        c12_fluxes[i].generateBinEdges(1E3, 2E6, 1000, EQUAL)

    py_printf('INFO', 'Simulating generational flux from H-1...')

    neutron = initializeNewNeutron()

    for i in range(num_neutrons):

        neutron._energy = max_energy				# initialize energy to 2 MeV

        for j in range(num_generations):
            h1.collideNeutron(neutron)
            h1_fluxes[j].tally(neutron)

    py_printf('INFO', 'Simulating generational flux from C-12...')

    for i in range(num_neutrons):

        neutron._energy = max_energy				# initialize energy to 2 MeV

        for j in range(num_generations):
            c12.collideNeutron(neutron)
            c12_fluxes[j].tally(neutron)			

    py_printf('INFO', 'Plotting the generational fluxes...')

    plotter.plotFluxes(h1_fluxes[1:], directory=output_dir, \
                            loglog=False, uselegend=False, \
                            filename='h1-generational-flux', \
                            title='H1 Generational Flux')

    plotter.plotFluxes(c12_fluxes[1:], directory=output_dir, \
                            loglog=False, uselegend=False, 
                            filename='c12-generational-flux', \
                            title='C12 Generational Flux')

    ###########################################################################
    #########################   Problems 3 and 4   ############################
    ###########################################################################

    num_generations = 50
    
    py_printf('NORMAL', 'Beginning Problems 3 and 4...')
    py_printf('NORMAL', '# generations = %d\t\t# neutrons / generation = %d', \
                                                num_generations, num_neutrons)
    py_printf('INFO', 'Initializing lethargy and energy flux tallies...')

    lethargy_flux = TallyFactory.createTally('lethargy flux', \
                                            c12, COLLISION_RATE)
    lethargy_flux.generateBinEdges(1E2, 2E6, 5000, LOGARITHMIC)

    energy_flux = TallyFactory.createTally('energy flux', c12, COLLISION_RATE)
    energy_flux.generateBinEdges(1E2, 2E6, 5000, EQUAL)

    neutron = initializeNewNeutron()

    py_printf('INFO', 'Simulating generational flux from C-12...')

    for i in range(num_neutrons):	
        neutron._energy = max_energy			# initialize energy to 2 MeV

        for j in range(num_generations):

            c12.collideNeutron(neutron)

            if (j > 0):
                lethargy_flux.tally(neutron)
                energy_flux.tally(neutron)
    
    py_printf('INFO', 'Plotting the lethargy and energy fluxes...')

    plotter.plotFlux(lethargy_flux, directory=output_dir, \
                        filename='lethargy-flux', \
                        title='C12 Lethargy Flux From Mono-energetic Source')
    plotter.plotFlux(energy_flux, directory=output_dir, \
                        filename='energy-flux', \
                        title='C12 Energy Flux From Mono-energetic Source')


    ###########################################################################
    #############################   Problem 5   ###############################
    ###########################################################################

    py_printf('NORMAL', 'Beginning Problem 5...')
    py_printf('NORMAL', '# generations = %d\t\t# neutrons / generation = %d', \
                                                num_generations, num_neutrons)

    num_generations = 15
    
    h1_material = Material('H-1')
    h1_material.setDensity(5., 'g/cc')
    h1_material.addIsotope(h1, 1.0)

    py_printf('INFO', 'Initializing lethargy flux tally...')

    lethargy_flux = TallyFactory.createTally('lethargy flux', \
                                            h1_material, FLUX)
    lethargy_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)

    neutron = initializeNewNeutron()
    fissioner = Fissioner()

    py_printf('INFO', 'Simulating H-1 flux for %d generations...', \
                                                        num_generations)

    for i in range(num_neutrons):	

        # Sample a fission energy from the Watt fission spectrum 
        neutron._energy = fissioner.emitNeutroneV()

        for j in range(num_generations):
            h1_material.collideNeutron(neutron)
            lethargy_flux.tally(neutron)
    
    py_printf('INFO', 'Plotting the lethargy flux...')

    plotter.plotFlux(lethargy_flux, directory=output_dir, \
                            filename='lethargy-flux-fission', \
                            title='H1 Energy Flux From Fission Spectrum')

    py_printf('NORMAL', 'Finished')


if __name__ == '__main__':

    main()  

