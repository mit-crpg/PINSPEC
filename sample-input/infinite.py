import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *

def main():
    
    # Set main simulation params
    num_batches = 10
    num_neutrons_per_batch = 10000
    num_threads = 4
    output_dir = 'Infinite'

    setLevel(INFO)
    
    # Call SLBW to create XS
    #py_printf('INFO', 'Creating SLBW xs')
    #T=300 #Temp in Kelvin of target nucleus
    #SLBW.replaceXS() #To clean previous changes to XS files
    #SLBW.SLBWXS('U-238',T,'capture') #To generate Doppler Broadened Res Cap
    #SLBW.SLBWXS('U-238',T,'scatter') #To generate Doppler Broadened Res Scat
    #SLBW.generatePotentialScattering('U-238') #To generate flat Res Scat XS
    #SLBW.compareXS('U-238', XStype='scatter', RI='no')
    #SLBW.compareXS('U-238', XStype='capture', RI='no')   

    # Define isotopes
    h1 = Isotope('H-1')
    b10 = Isotope('B-10')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    zr90 = Isotope('Zr-90')    

    # Define materials
    mix = Material('Fuel Moderator Mix')
    mix.setDensity(5., 'g/cc')
    #mix.addIsotope(b10, .0000001)
    #mix.addIsotope(o16, 1.0)
    mix.addIsotope(h1, 1.0)
    mix.addIsotope(u238, 0.0000010)
    #mix.addIsotope(u235, .0025)
    #mix.addIsotope(zr90, .16)

    py_printf('INFO', 'Added isotopes')
    
    # Define regions
    region_mix = Region('infinite medium', INFINITE)
    region_mix.setMaterial(mix)

    py_printf('INFO', 'Made mixture region')

    # Create a tally for the flux
 
    # Create a tallies to compute the the RI
    abs_rate = TallyFactory.createTally('absorption rate', \
												u238, ABSORPTION_RATE)
    flux_RI = TallyFactory.createTally('flux RI', mix, FLUX)
    RI_bin_edges = numpy.array([0.01, 0.1, 1.0, 6.0, 10.0, 25.0, \
												50.0, 100.0, 1000.0])
    abs_rate.setBinEdges(RI_bin_edges)
    flux_RI.setBinEdges(RI_bin_edges)

    # Set a precision trigger: tells simulation to run until maximum relative
    # error is less than the trigger value (4E-3)
    abs_rate.setPrecisionTrigger(RELATIVE_ERROR, 2E-2)

    # Define geometry
    geometry = Geometry()
    geometry.setSpatialType(INFINITE_HOMOGENEOUS)
    geometry.addRegion(region_mix)
    flux = TallyFactory.createTally('total flux', region_mix, FLUX)
    flux.generateBinEdges(1E-2, 1E7, 10000, LOGARITHMIC) 

	# Register the tallies
    TallyBank.registerTally(flux, geometry)
    TallyBank.registerTally(abs_rate, mix)
    TallyBank.registerTally(flux_RI, mix)

    geometry.setNumBatches(num_batches)
    geometry.setNeutronsPerBatch(num_neutrons_per_batch)
    geometry.setNumThreads(num_threads)

    py_printf('INFO', 'Made geometry')

    # Run Monte Carlo simulation
    geometry.runMonteCarloSimulation();


    # Plot the microscopic cross sections for each isotope
    plotter.plotMicroXS(u235, ['capture', 'elastic', 'fission'], output_dir)
    plotter.plotMicroXS(u238, ['capture', 'elastic', 'fission'], output_dir)
    plotter.plotMicroXS(h1, ['capture', 'elastic', 'absorption'], output_dir)
    plotter.plotMicroXS(o16, ['capture', 'elastic', 'absorption'], output_dir)
    plotter.plotMacroXS(mix, ['total', 'capture', 'elastic', \
                                        'fission', 'absorption'], output_dir)

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
    TallyBank.outputBatchStatistics(output_dir, 'test')

    # plot the fission spectrum the CDF
    plotter.plotFissionSpectrum(output_dir)

    #Plot the thermal scattering kernel PDFs and CDFs
    plotter.plotThermalScattering(h1, output_dir)

if __name__ == '__main__':

    main()  

