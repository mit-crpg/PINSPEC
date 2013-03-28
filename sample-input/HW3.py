# This input file corresponds to 22.211 HW3, where we have: 
# slowing down MC, U238 absorption resonance xs (and 0.1 above the 
# energy range). Scattering xs is asymptotic down-scattering. 
import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.log as log

def main():

    # Set main simulation params
    num_batches            = 10
    num_neutrons_per_batch = 10000
    num_threads            = 4
    output_dir             = 'HW3'
    
    # set logging level
    log.setLevel('INFO')
    
    # make geometry
    geometry = Geometry()

    # Call SLBW to create XS
    filename = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
    Temp = 300 #Temp in Kelvin of target nucleus
    SLBW.SLBWXS('U-238',Temp,'capture') #To generate Doppler Broadened Res Cap
    SLBW.compareXS('U-238', XStype='capture', RI='no')   
    
    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
 
    # Define materials
    mix = Material()
    mix.setMaterialName('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc') # we do not have a density from last year.
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(u238, 0.40)
    mix.addIsotope(u235, .02) # we do not have U235 from last year. 

    log.py_printf(INFO, 'Added isotopes')

    # Define regions
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    log.py_printf(INFO, 'Made mixture region')
    
    # Define geometry
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    # Create a tally for the flux
    flux = Tally('total flux', GEOMETRY, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC)
    geometry.addTally(flux)
    
    # Create a tally for the absorption rate
    abs_rate = Tally('absorption rate', region_mix, ABSORPTION_RATE)
    bin_edges = numpy.array([0.1, 1., 5., 10., 100., 1000.])
    abs_rate.setBinEdges(bin_edges)
    geometry.addTally(abs_rate)
    abs_rate.setPrecisionTrigger(RELATIVE_ERROR, 3E-3)
    
    log_printf(INFO, 'Made geometry')
   
    # Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();
    
    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics(output_dir, 'test')

    # Plotting  
    plotter.plotFlux(flux, output_dir)
    plotter.plotMicroXS(u235, ['capture', 'absorption'])
    plotter.plotMicroXS(u238, ['capture', 'fission', 'elastic', \
                                   'total', 'absorption'])
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'])
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'])
    plotter.plotMacroXS(mix, ['capture', 'elastic', 'fission', \
                                  'absorption', 'total'])
    
    # plot the fission spectrum the CDF
    plotter.plotFissionSpectrum()

    #Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScattering(h1, output_dir)
    
    

if __name__ == '__main__':

    main()  

