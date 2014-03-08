import numpy
from pinspec import *
import pinspec.plotter as plotter
from pinspec.log import *


###############################################################################
##    This python script uses PINSPEC to generate results for homework 1 of
##    the 22.211 Introduction to Reactor Physics course taught at MIT by
##    Prof. Kord Smith in spring 2012.
###############################################################################


###############################################################################
###########################  Main Simulation Paramters  #######################
###############################################################################

# Set main simulation params
num_neutrons = 100000
setOutputDirectory('generational');
py_setlevel('INFO')

py_printf('TITLE', 'Starting a simulation to compute generational fluxes')


###############################################################################
##############################  Create Isotopes ###############################
###############################################################################

py_printf('NORMAL', 'Initializing isotopes...')
h1 = Isotope('H-1')
c12 = Isotope('C-12')

# Zero out the capture xs
energy = numpy.array([1E-7, 2e7])	        # energy bounds
xs = numpy.array([0.0])				# one group xs
h1.setMultigroupCaptureXS(energy, xs)   
c12.setMultigroupCaptureXS(energy, xs)

py_printf('INFO', 'Plotting microscopic cross-sections...')
plotter.plotMicroXS(h1, ['capture', 'elastic', 'fission'])
plotter.plotMicroXS(c12, ['capture', 'elastic', 'fission'])

# Turn off thermal scattering
h1.neglectThermalScattering()
c12.neglectThermalScattering()


###############################################################################
###########################   Problems 1 and 2   ##############################
###############################################################################

num_generations = 15
max_energy = 2E6

py_printf('HEADER', 'Problems 1 and 2')
py_printf('INFO', 'Initializing generational flux tallies...')

# Create flux tallies for each generation
h1_fluxes = []
c12_fluxes = []

for i in range(num_generations):

    h1_fluxes.append(TallyFactory.createTally(h1, COLLISION_RATE))
    h1_fluxes[i].generateBinEdges(0., 2E6, 1000, EQUAL)

    c12_fluxes.append(TallyFactory.createTally(c12, COLLISION_RATE))
    c12_fluxes[i].generateBinEdges(1E3, 2E6, 1000, EQUAL)

py_printf('INFO', 'Simulating generational flux from H-1...')
py_printf('INFO', '# neutrons = %d\t\t# generations = %d', 
                                            num_neutrons, num_generations)

neutron = createNewNeutron()

for i in range(num_neutrons):

    neutron._energy = max_energy	     # initialize energy to 2 MeV

    for j in range(num_generations):
        h1.collideNeutron(neutron)
        h1_fluxes[j].tally(neutron)

py_printf('INFO', 'Simulating generational flux from C-12...')

for i in range(num_neutrons):

    neutron._energy = max_energy	     # initialize energy to 2 MeV

    for j in range(num_generations):
        c12.collideNeutron(neutron)
        c12_fluxes[j].tally(neutron)			

py_printf('INFO', 'Plotting the generational fluxes...')

plotter.plotFlux(h1_fluxes[1:], loglog=False, uselegend=False, \
                  filename='h1-generational-flux', title='H1 Generational Flux')

plotter.plotFlux(c12_fluxes[1:], loglog=False, uselegend=False, 
                filename='c12-generational-flux', title='C12 Generational Flux')


###############################################################################
###########################   Problems 3 and 4   ##############################
###############################################################################

num_generations = 50
    
py_printf('HEADER', 'Problems 3 and 4')
py_printf('INFO', 'Initializing lethargy and energy flux tallies...')

lethargy_flux = TallyFactory.createTally(c12, COLLISION_RATE)
lethargy_flux.generateBinEdges(1E2, 2E6, 5000, LOGARITHMIC)

energy_flux = TallyFactory.createTally(c12, COLLISION_RATE)
energy_flux.generateBinEdges(1E2, 2E6, 5000, EQUAL)

neutron = createNewNeutron()

py_printf('INFO', 'Simulating generational flux from C-12...')
py_printf('INFO', '# neutrons = %d\t\t# generations = %d', 
                                                num_neutrons, num_generations)
for i in range(num_neutrons):	
    neutron._energy = max_energy	     # initialize energy to 2 MeV

    for j in range(num_generations):
        c12.collideNeutron(neutron)

        if (j > 0):
            lethargy_flux.tally(neutron)
            energy_flux.tally(neutron)
    
py_printf('INFO', 'Plotting the lethargy and energy fluxes...')

plotter.plotFlux(lethargy_flux, filename='lethargy-flux', \
                        title='C12 Lethargy Flux From Mono-energetic Source')
plotter.plotFlux(energy_flux, filename='energy-flux', \
                        title='C12 Energy Flux From Mono-energetic Source')


###############################################################################
###############################   Problem 5   #################################
###############################################################################

py_printf('HEADER', 'Problem 5')

num_generations = 15
    
h1_material = Material('H-1')
h1_material.setDensity(0.07778, 'g/cc')
h1_material.addIsotope(h1, 2.0)

py_printf('INFO', 'Initializing lethargy flux tally...')

lethargy_flux = TallyFactory.createTally(h1_material, FLUX)
lethargy_flux.generateBinEdges(1E-2, 1E7, 5000, LOGARITHMIC)

neutron = createNewNeutron()
fissioner = Fissioner()

py_printf('INFO', 'Simulating H-1 flux for %d generations...', \
                                                        num_generations)
py_printf('INFO', '# neutrons = %d\t\t# generations = %d', 
                                                num_neutrons, num_generations)

for i in range(num_neutrons):	

    # Sample a fission energy from the Watt fission spectrum 
    neutron._energy = fissioner.emitNeutroneV()

    for j in range(num_generations):
        h1_material.collideNeutron(neutron)
        lethargy_flux.tally(neutron)
    
py_printf('INFO', 'Plotting the lethargy flux...')

plotter.plotFlux(lethargy_flux, filename='lethargy-flux-fission', \
                            title='H1 Energy Flux From Fission Spectrum')

py_printf('TITLE', 'Finished')
