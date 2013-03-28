import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
import pinspec.log as log

def main():
    
    # Set main simulation params
    num_batches            = 10
    num_neutrons_per_batch = 10000
    num_threads            = 4
    output_dir             = 'Infinite'

    # set logging level
    log.setLevel('INFO')
    
    # make geometry
    geometry = Geometry()
    
    log.py_printf('INFO', 'Creating SLBW xs')    
    
    # Call SLBW to create XS    
    Temp = 300 #Temp in Kelvin of target nucleus
    SLBW.replaceXS() #To clean previous changes to XS files
    SLBW.SLBWXS('U-238',Temp,'capture') #To generate Doppler Broadened Res Cap
    #SLBW.SLBWXS('U-238',Temp,'scatter') #To generate Doppler Broadened Res Scat
    #SLBW.generatePotentialScattering('U-238') #To generate flat Res Scat XS
    #SLBW.compareXS('U-238', XStype='scatter', RI='no')
    SLBW.compareXS('U-238', XStype='capture', RI='no')   

    # Define isotopes
    h1   = Isotope('H-1')
    b10  = Isotope('B-10')
    o16  = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    zr90 = Isotope('Zr-90')    

    # Plot the microscopic cross sections for each isotope
    plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission'], output_dir)
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'], output_dir)
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'], output_dir)
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'], output_dir)
    plotter.plotMacroXS(mix, ['total', 'capture', 'elastic',
                              'fission', 'absorption'], output_dir)

    # Define materials
    mix = Material()
    mix.setMaterialName('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc')
    mix.addIsotope(o16, 1.0)
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(u238, 0.40)
    mix.addIsotope(u235, .02)
    mix.addIsotope(zr90, .16)

    log.py_printf('INFO', 'Added isotopes')
    
    # Define regions
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    log.py_printf('INFO', 'Made mixture region')

    # Define geometry
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)
    geometry.addRegion(region_mix)

    # create and define tallies
    abs_rate  = Tally('absorption rate', u238, ABSORPTION_RATE)
    flux_RI   = Tally('flux RI', region_mix, FLUX)
    flux      = Tally('total flux', geometry, FLUX)
    bin_edges = numpy.array([0.1, 1., 6., 10., 25., 50., 100.])
    abs_rate.setBinEdges(bin_edges)
    flux_RI.setBinEdges(bin_edges)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC) 
    
    # add tallies to geometry
    geometry.addTally(flux)
    geometry.addTally(abs_rate)
    geometry.addTally(flux_RI)
    
    log.py_printf('INFO', 'Geometry setup complete')

    # Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();

    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics(output_dir, 'test')

    # Plot fluxes
    plotter.plotFlux(flux, output_dir)
    
    # Compute the resonance integrals
    RI = process.RI(flux_RI, abs_rate)
    RI.printRI()
    plotter.plotRI(RI, output_dir)
    
    # Compute the group cross sections
    groupXS = process.groupXS(flux_RI, abs_rate)
    groupXS.printGroupXS()
    plotter.plotGroupXS(groupXS, output_dir)

    # Dump batch statistics to output files to some new directory
    geometry.outputBatchStatistics(output_dir, 'test')

    # Plot the fission spectrum the CDF
    plotter.plotFissionSpectrum(output_dir)

    # Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScattering(h1, output_dir)

if __name__ == '__main__':

    main()  

