import numpy
from pinspec import *
import pinspec.plotter as plotter
from pinspec.log import *


###############################################################################
##    This python script uses PINSPEC to generate a flux for an arbitrary
##    infinite medium, containing the predominant isotopes found in a PWR.
###############################################################################

setOutputDirectory('infinite')

py_setlevel('INFO')
py_printf('TITLE', 'Simulating an infinite medium homogenized pin-cell')
    

###############################################################################
##############################  Create Isotopes ###############################
###############################################################################

py_printf('NORMAL', 'Initializing isotopes...')
h1 = Isotope('H-1')
b10 = Isotope('B-10')
o16 = Isotope('O-16')
u235 = Isotope('U-235')
u238 = Isotope('U-238')
zr90 = Isotope('Zr-90')
    

###############################################################################
##############################  Create Materials ##############################
###############################################################################

py_printf('NORMAL', 'Initializing infinite medium material...')
mix = Material('Fuel Moderator Mix')
mix.setDensity(5., 'g/cc')
mix.addIsotope(b10, .0000001)
mix.addIsotope(o16, 1.0)
mix.addIsotope(h1, 1.0)
mix.addIsotope(u238, 0.01)
mix.addIsotope(u235, .0025)
mix.addIsotope(zr90, .16)


###############################################################################
###############################  Create Regions ###############################
###############################################################################

py_printf('NORMAL', 'Initializing infnite medium region...')
region_mix = RegionFactory.createRegion(INFINITE_MEDIUM)
region_mix.setMaterial(mix)


###############################################################################
###############################  Create Geometry ##############################
###############################################################################
    
py_printf('NORMAL', 'Initializing the geometry...')
geometry = Geometry(INFINITE_HOMOGENEOUS)
geometry.addRegion(region_mix)


###############################################################################
################################  Create Tallies ##############################
###############################################################################
    
py_printf('NORMAL', 'Initializing flux tally...')
flux = TallyFactory.createTally(region_mix, FLUX)
flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC)
TallyBank.registerTally(flux)
   

###############################################################################
#########################  Run Monte Carlo Simulation #########################
###############################################################################

# Run Monte Carlo simulation
geometry.runMonteCarloSimulation()


###############################################################################
############################  Process Output Data #############################
###############################################################################

py_printf('INFO', 'Plotting flux...')
plotter.plotFlux(flux, uselegend=False, filename='flux', \
                     title='Infinite medium flux')
    
py_printf('INFO', 'Writing tally batch statistics to output file...')
TallyBank.outputBatchStatistics()

    
py_printf('TITLE', 'Finished')
