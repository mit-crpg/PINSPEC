import numpy
from pinspec import *
from pinspec.log import *


###############################################################################
##    This python script uses PINSPEC to generate results for an
##    outscattering rate tally.
###############################################################################



###############################################################################
###########################  Main Simulation Paramters  #######################
###############################################################################

set_output_directory('outscatter')

py_set_log_level('INFO')

py_printf('TITLE', 'Simulation to tally outscattering')


###############################################################################
##############################  Create Isotopes ###############################
###############################################################################

py_printf('NORMAL', 'Initializing isotopes...')
h1 = Isotope('H-1')
o16 = Isotope('O-16')
u235 = Isotope('U-235')
u238 = Isotope('U-238')
    

###############################################################################
##############################  Create Materials ##############################
###############################################################################

py_printf('NORMAL', 'Initializing fuel-moderator mix material...')
mix = Material('Fuel Moderator Mix')
mix.setDensity(5., 'g/cc')
mix.addIsotope(o16, 1.0)
mix.addIsotope(h1, 1.0)
mix.addIsotope(u238, 0.01)
mix.addIsotope(u235, .0025)
    

###############################################################################
###############################  Create Regions ###############################
###############################################################################

py_printf('NORMAL', 'Initializing fuel-moderator mix region...')
region_mix = InfiniteMediumRegion('infinite medium')
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

py_printf('NORMAL', 'Initializing tally...')
outscatter_rate = TallyFactory.createTally(region_mix, OUTSCATTER_RATE, \
                                               'outscatter rate')
outscatter_rate.generateBinEdges(1E-2, 1E7, 5, LOGARITHMIC)
TallyBank.registerTally(outscatter_rate)


###############################################################################
#########################  Run Monte Carlo Simulation #########################
###############################################################################

# Run Monte Carlo simulation
geometry.runMonteCarloSimulation()


###############################################################################
############################  Process Output Data #############################
###############################################################################

py_printf('INFO', 'Writing tally batch statistics to output file...')
TallyBank.outputBatchStatistics()
outscatter_rate.printTallies()
    
py_printf('TITLE', 'Finished')
