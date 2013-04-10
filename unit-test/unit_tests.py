import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *
import unittest
import os


# 1) Unit test important C++ and python functions
    
class TestSimulationParameters(unittest.TestCase):
    
    
    def testLogLevel(self):
        py_printf('UNITTEST', 'Testing Simulation Parameters py_setlevel')
        py_setlevel('UNITTEST')
        log_val = assignValue('UNITTEST')
        self.assertEqual(log_val, get_loglevel())
    
    
    def testSetOutputDirectory(self):
        py_printf('UNITTEST', 'Testing Simulation Parameters setOutputDirectory')
        setOutputDirectory('infinite')
        self.assertEqual(getOutputDirectory(),'infinite')


class TestXSReader(unittest.TestCase):
    
    
    # check to make sure restoreXSLibrary function is able to execute function to copy backup xs files
    def testRestoreXSLibrary(self):
        py_printf('UNITTEST', 'Testing xsreader restoreXSLibrary')
        ret = restoreXSLibrary()
        self.assertEqual(ret, 0)

    
    # check to make sure all xs files return at least 1 xs data point
    def testGetNumCrossSectionDataPoints(self):
        py_printf('UNITTEST', 'Testing xsreader GetNumCrossSectionDataPoints')
        pkg_path = getXSLibDirectory()
        xs_paths = os.listdir(pkg_path)
        
        for xs_path in xs_paths:
            if xs_path[-4:] == '.txt':
                xs_full_path = os.path.join(xs_lib_path, xs_path)
                num_points = getNumCrossSectionDataPoints(xs_full_path)
                self.assertGreater(num_points, 0)


class TestIsotope(unittest.TestCase):
    
    
    # Set main simulation params
    def setUp(self):
        self.num_batches = 5
        self.num_neutrons_per_batch = 10000
        self.num_threads = 4

    
    # Test Isotope constructor
    def testIsotope(self):
        py_printf('UNITTEST', 'Testing Isotope Constructor')
        try:
            h1 = Isotope('H-1')
        except:
            self.fail('Isotope constructor function failed')

                
    # Test IsotopeName function
    def testIsotopeName(self):
        py_printf('UNITTEST', 'Testing Isotope IsotopeName')
        h1 = Isotope('H-1')
        self.assertEqual(h1.getIsotopeName(), 'H-1')

        
    # Test Isotope constructor to see if Alpha is set to the correct value
    def testAlpha(self):
        py_printf('UNITTEST', 'Testing Isotope set Alpha')
        o16 = Isotope('O-16')
        self.assertAlmostEqual(o16.getAlpha(), (15.0/17)**2)

    
    # Test that a fissionable isotope is set to fissionable
    def testFissionable(self):
        py_printf('UNITTEST', 'Testing Isotope isFissionable')
        u235 = Isotope('U-235')
        self.assertTrue(u235.isFissionable())

    
    # Test setTemperature function
    def testTemperature(self):
        py_printf('UNITTEST', 'Testing Isotope setTemperature')
        o16 = Isotope('O-16')
        o16.setTemperature(500)
        self.assertEqual(o16.getTemperature(), 500)


    # Test Isotope constructor to see if A is set to the correct value
    def testA(self):
        py_printf('UNITTEST', 'Testing Isotope set A')
        o16 = Isotope('O-16')
        self.assertEqual(o16.getA(), 16)


    # Test Isotope setAO function
    def testAO(self):
        py_printf('UNITTEST', 'Testing Isotope setAO')
        o16 = Isotope('O-16')
        o16.setAO(1.1)
        self.assertAlmostEqual(o16.getAO(), 1.1)

    
    # Test Isotope setN function
    def testN(self):
        py_printf('UNITTEST', 'Testing Isotope setN')
        o16 = Isotope('O-16')
        o16.setN(1.1)
        self.assertAlmostEqual(o16.getN(), 1.1)

    
    # Test Isotope getMuAverage function
    def testMuAverage(self):
        py_printf('UNITTEST', 'Testing Isotope getMuAverage')
        o16 = Isotope('O-16')
        self.assertAlmostEqual(o16.getMuAverage(), 2.0/(3.0*16))
    
    
    # Test Isotope getNumXSEnergies function
    def testGetNumXSEnergies(self):
        py_printf('UNITTEST', 'Testing Isotope getNumXSEnergies')
        u235 = Isotope('U-235')
        self.assertLess(0, u235.getNumXSEnergies('elastic'))
    
    
    # Test Isotope retrieveXSEnergies function
    def testRetrieveXSEnergies(self):
        py_printf('UNITTEST', 'Testing Isotope retrieveXSEnergies')
        u235 = Isotope('U-235')
        num_energies = u235.getNumXSEnergies('elastic')
        energies = u235.retrieveXSEnergies(num_energies, 'elastic')
        self.assertEqual(len(energies), num_energies)


    # Test Isotope retrieveXS function
    def testRetrieveXS(self):
        py_printf('UNITTEST', 'Testing Isotope retrieveXS')
        u235 = Isotope('U-235')
        num_energies = u235.getNumXSEnergies('elastic')
        xs = u235.retrieveXS(num_energies, 'elastic')
        self.assertEqual(len(xs), num_energies)

    
    # Test Isotope get xs functions with Energy input
    def testGetXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Isotope get xs functions by Energy')
        u235 = Isotope('U-235')
        self.assertLessEqual(abs((u235.getTotalXS(1.0e3) + u235.getFissionXS(1.0e3) + u235.getTransportXS(1.0e3)) / (u235.getElasticXS(1.0e3) + u235.getAbsorptionXS(1.0e3) + u235.getFissionXS(1.0e3) + + u235.getTransportXS(1.0e3)) - 1), .001)


    # Test Isotope get xs functions with Index input
    def testGetXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Isotope get xs functions by Index')
        u235 = Isotope('U-235')
        self.assertLessEqual(abs((u235.getTotalXS(1e3) + u235.getFissionXS(1e3) + u235.getTransportXS(1e3)) / (u235.getElasticXS(1e3) + u235.getAbsorptionXS(1e3) + u235.getFissionXS(1e3)+ u235.getTransportXS(1e3)) - 1), .001)
    
    
    # Test Isotope setElasticXS function
    def testsetElasticXS(self):
        py_printf('UNITTEST', 'Testing Isotope setElasticXS')
        u235 = Isotope('U-235')
        energies = numpy.array([1E-7, 2E7])
        xs = numpy.array([1, 2])
        
        try:
            u235.setElasticXS(energies, xs)
        except:
            self.fail('Unable to set elastic XS')


    # Test Isotope setMultigroupElasticXS function
    def testSetMultigroupElasticXS(self):
        py_printf('UNITTEST', 'Testing Isotope setMultigroupElasticXS')
        u235 = Isotope('U-235')
        energies = numpy.array([1E-7, 1E0, 2E7])
        xs = numpy.array([1, 2])
        
        try:
            u235.setMultigroupElasticXS(energies, xs)
        except:
            self.fail('Unable to set elastic multigroup XS')


    # Test Isotope setFissionXS function
    def testsetElasticXS(self):
        py_printf('UNITTEST', 'Testing Isotope setFissionXS')
        u235 = Isotope('U-235')
        energies = numpy.array([1E-7, 2E7])
        xs = numpy.array([1, 2])
        
        try:
            u235.setFissionXS(energies, xs)
        except:
            self.fail('Unable to set fission XS')
    
    
    # Test Isotope setMultigroupFissionXS function
    def testSetMultigroupFissionXS(self):
        py_printf('UNITTEST', 'Testing Isotope setMultigroupFissionXS')
        u235 = Isotope('U-235')
        energies = numpy.array([1E-7, 1E0, 2E7])
        xs = numpy.array([1, 2])
        
        try:
            u235.setMultigroupFissionXS(energies, xs)
        except:
            self.fail('Unable to set fission XS')



    # Test Isotope setCaptureXS function
    def testsetCaptureXS(self):
        py_printf('UNITTEST', 'Testing Isotope setCaptureXS')
        u235 = Isotope('U-235')
        energies = numpy.array([1E-7, 2E7])
        xs = numpy.array([1, 2])
        
        try:
            u235.setCaptureXS(energies, xs)
        except:
            self.fail('Unable to set capture XS')


    # Test Isotope setMultigroupCaptureXS function
    def testSetMultigroupCaptureXS(self):
        py_printf('UNITTEST', 'Testing Isotope setMultigroupCaptureXS')
        u235 = Isotope('U-235')
        energies = numpy.array([1E-7, 1E0, 2E7])
        xs = numpy.array([1, 2])
        
        try:
            u235.setMultigroupCaptureXS(energies, xs)
        except:
            self.fail('Unable to set capture XS')


    # Test Isotope neglectThermalScattering function
    def testNeglectThermalScattering(self):
        py_printf('UNITTEST', 'Testing Isotope neglectThermalScattering')
        u235 = Isotope('U-235')
        u235.neglectThermalScattering()
        self.assertFalse(u235.usesThermalScattering())


    # Test Isotope useThermalScattering function
    def testUseThermalScattering(self):
        py_printf('UNITTEST', 'Testing Isotope useThermalScattering')
        u235 = Isotope('U-235')
        u235.useThermalScattering()
        self.assertTrue(u235.usesThermalScattering())


    # Test Isotope isRescaled function
    def testIsRescaled(self):
        py_printf('UNITTEST', 'Testing Isotope isRescaled')
        u235 = Isotope('U-235')
        self.assertTrue(u235.isRescaled())


    # Test Isotope isRescaled function
    def testSetA(self):
        py_printf('UNITTEST', 'Testing Isotope setA')
        u235 = Isotope('U-235')
        u235.setA(10)
        self.assertEqual(u235.getA(),10)


    # Test Isotope isRescaled function
    def testSetA(self):
        py_printf('UNITTEST', 'Testing Isotope setA')
        u235 = Isotope('U-235')
        u235.setA(10)
        self.assertEqual(u235.getA(),10)


    # Test Isotope clone function
    def testClone(self):
        py_printf('UNITTEST', 'Testing Isotope clone')
        u235 = Isotope('U-235')
        u235_clone = u235.clone()
        self.assertEqual(u235.getA(),u235_clone.getA())


    # Test Isotope getThermalScatteringEnergy function
    def testGetThermalScatteringEnergy(self):
        py_printf('UNITTEST', 'Testing Isotope getThermalScatteringEnergy')
        u235 = Isotope('U-235')
        self.assertGreater(u235.getThermalScatteringEnergy(1.0),0.0)


    # Test Isotope getNumThermalCDFs function
    def testgetNumThermalCDFs(self):
        py_printf('UNITTEST', 'Testing Isotope getNumThermalCDFs')
        u235 = Isotope('U-235')
        self.assertGreater(u235.getNumThermalCDFs(),0)


    # Test Isotope getNumThermalCDFBins function
    def testgetNumThermalCDFBins(self):
        py_printf('UNITTEST', 'Testing Isotope getNumThermalCDFBins')
        u235 = Isotope('U-235')
        self.assertGreater(u235.getNumThermalCDFBins(),0)


    # Test Isotope retrieveEprimeToE function
    def testRetrieveEprimeToE(self):
        py_printf('UNITTEST', 'Testing Isotope retrieveEprimeToE')
        h1 = Isotope('H-1')
        num_bins = h1.getNumThermalCDFBins()
        num_cdfs = h1.getNumThermalCDFs()
        Eprime_to_E = h1.retrieveEprimeToE(num_bins)
        self.assertEqual(num_bins,len(Eprime_to_E))


    # Test Isotope retrieveEtokT function
    def testRetrieveEtokT(self):
        py_printf('UNITTEST', 'Testing Isotope retrieveEtokT')
        h1 = Isotope('H-1')
        num_bins = h1.getNumThermalCDFBins()
        num_cdfs = h1.getNumThermalCDFs()
        E_to_kT = h1.retrieveEtokT(num_cdfs)
        self.assertEqual(num_cdfs,len(E_to_kT))

        
    # Test Isotope retrieveThermalCDFs function
    def testRetrieveThermalCDFs(self):
        py_printf('UNITTEST', 'Testing Isotope retrieveThermalCDFs')
        h1 = Isotope('H-1')
        num_bins = h1.getNumThermalCDFBins()
        num_cdfs = h1.getNumThermalCDFs()
        cdfs = h1.retrieveThermalCDFs(num_cdfs*num_bins)
        self.assertEqual(num_cdfs*num_bins,len(cdfs))
        

    # Test Isotope retrieveThermalDistributions function
    def testRetrieveThermalDistributions(self):
        py_printf('UNITTEST', 'Testing Isotope retrieveEprimeToE')
        h1 = Isotope('H-1')
        num_bins = h1.getNumThermalCDFBins()
        num_cdfs = h1.getNumThermalCDFs()
        dist = h1.retrieveThermalDistributions(num_cdfs*num_bins)
        self.assertEqual(num_cdfs*num_bins,len(dist))

        
    # Test Isotope getDistanceTraveled function
    def testGetDistanceTraveled(self):
        py_printf('UNITTEST', 'Testing Isotope getDistanceTraveled')
        h1 = Isotope('H-1')
        neutron = initializeNewNeutron()
        neutron._energy = 1.0
        dist = h1.getDistanceTraveled(neutron)
        self.assertGreater(dist, 0.0)


    # Test Isotope collideNuetron function
    def testCollideNeutron(self):
        py_printf('UNITTEST', 'Testing Isotope collideNeutron')
        h1 = Isotope('H-1')
        neutron = initializeNewNeutron()
        neutron._energy = 1.0
        dist = h1.collideNeutron(neutron)
        self.assertEqual(neutron._old_energy, 1.0)


class TestMaterial(unittest.TestCase):

    def setUp(self):

        # Set main simulation params
        self.num_batches = 5
        self.num_neutrons_per_batch = 10000
        self.num_threads = 4
    
        # Define isotopes
        self.h1 = Isotope('H-1')
        self.o16 = Isotope('O-16')


    # Test Material getMaterialName function
    def testGetMaterialName(self):
        py_printf('UNITTEST', 'Testing Material getMaterialName')
        mod = Material('mod')
        self.assertEqual(mod.getMaterialName(), 'mod')


    # Test Material getMaterialNumberDensity function
    def testGetMaterialNumberDensity(self):
        py_printf('UNITTEST', 'Testing Material getMaterialNumberDensity')
        mod = Material('mod')
        mod.setDensity(1., 'at/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertAlmostEqual(mod.getMaterialNumberDensity(), 1.0)

    
    # Test Material getIsotopeNumDensity function
    def testGetIsotopeNumDensity(self):
        py_printf('UNITTEST', 'Testing Material getIsotopeNumDensity')
        mod = Material('mod')
        mod.setDensity(1., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertAlmostEqual(mod.getIsotopeNumDensity(self.h1) / 6.02299993035875E23 , 1.0)
    
    
    # Test Material containsIsotope function
    def testContainsIsotope(self):
        py_printf('UNITTEST', 'Testing Material containsIsotope')
        mod = Material('mod')
        mod.setDensity(1., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertTrue(mod.containsIsotope(self.h1))
    
    
    # Test Material getVolume function
    def testGetVolume(self):
        py_printf('UNITTEST', 'Testing Material getVolume')
        mod = Material('mod')
        mod.incrementVolume(5.)
        self.assertEqual(mod.getVolume(), 5.)


    # Test Material getBucklingSquared function
    def testGetBucklingSquared(self):
        py_printf('UNITTEST', 'Testing Material getBucklingSquared')
        mod = Material('mod')
        mod.setBucklingSquared(5.)
        self.assertEqual(mod.getBucklingSquared(), 5.)

    
    # Test Material getNumXSEnergies function
    def testGetNumXSEnergies(self):
        py_printf('UNITTEST', 'Testing Material getNumXSEnergies')
        mod = Material('mod')
        mod.setDensity(1., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertGreater(mod.getNumXSEnergies('elastic'),0)
    

    # Test Material getDensity function
    def testGetDensity(self):
        py_printf('UNITTEST', 'Testing Material getDensity')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertEqual(mod.getDensity(), 5.)

    
    # Test Material containsIsotope function
    def testContainsIsotope(self):
        py_printf('UNITTEST', 'Testing Material containsIsotope')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertTrue(mod.containsIsotope(self.h1))


    # Test Material getXSbyEnergy function
    def testGetXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getXSbyEnergy')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertLessEqual(abs(mod.getTotalMacroXS(1.0e3) / (mod.getElasticMacroXS(1.0e3) + mod.getAbsorptionMacroXS(1.0e3)) - 1), .001)


    # Test Material getXSbyIndex function
    def testGetXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Material getXSbyIndex')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertLessEqual(abs(mod.getTotalMacroXS(1e3) / (mod.getElasticMacroXS(1e3) + mod.getAbsorptionMacroXS(1e3)) - 1), .001)

    # Test Material retrieveXSEnergies function
    def testRetrieveXSEnergies(self):
        py_printf('UNITTEST', 'Testing Material retrieveXSEnergies')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        num_energies = 0
        num_energies = mod.getNumXSEnergies('elastic')
        energies = mod.retrieveXSEnergies(num_energies, 'elastic')
        self.assertEqual(len(energies), num_energies)


    # Test Material retrieveXS function
    def testRetrieveXS(self):
        py_printf('UNITTEST', 'Testing Material retrieveXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        num_energies = 0
        num_energies = mod.getNumXSEnergies('elastic')
        xs = mod.retrieveXS(num_energies, 'elastic')
        self.assertEqual(len(xs), num_energies)


    # Test Material getTotalMacroXS function
    def testGetTotalMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getTotalMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getTotalMacroXS(1.0), 0.0)


    # Test Material getTotalMacroXS function
    def testGetTotalMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Material getTotalMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getTotalMacroXS(100), 0.0)


    # Test Material getTotalMacroXS function
    def testGetTotalMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getTotalMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getTotalMacroXS(1.0), 0.0)


    # Test Material getElasicMacroXS function
    def testGetElasticMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Material getElasticMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getElasticMacroXS(100), 0.0)
    
    
    # Test Material getElasticMacroXS function
    def testGetTotalMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getElasticMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getElasticMacroXS(1.0), 0.0)

    
    # Test Material getAbsorptionMacroXS function
    def testGetAbsorptionMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Material getAbsorptionMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getAbsorptionMacroXS(100), 0.0)
    
    
    # Test Material getAbsorptionMacroXS function
    def testGetAbsorptionMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getAbsorptionMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getAbsorptionMacroXS(1.0), 0.0)


    # Test Material getCaptureMacroXS function
    def testGetCaptureMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Material getCaptureMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getCaptureMacroXS(100), 0.0)
    
    
    # Test Material getTotalMacroXS function
    def testGetCaptureMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getCaptureMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getCaptureMacroXS(1.0), 0.0)


    # Test Material getTransportMacroXS function
    def testGetTransportMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Material getTransportMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getTransportMacroXS(100), 0.0)
    
    
    # Test Material getTotalMacroXS function
    def testGetTransportMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Material getTransportMacroXS')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod.addIsotope(self.o16, 0.5)
        self.assertGreater(mod.getTransportMacroXS(1.0), 0.0)


    # Test Material addIsotope function
    def testAddIsotope(self):
        py_printf('UNITTEST', 'Testing Material addIsotope')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertTrue(mod.containsIsotope(self.h1))


    # Test Material addIsotope function
    def testUnAddIsotope(self):
        py_printf('UNITTEST', 'Testing Material un-addIsotope')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        self.assertFalse(mod.containsIsotope(self.o16))


    # Test Material sampleIsotope function
    def testSampleIsotope(self):
        py_printf('UNITTEST', 'Testing Material sampleIsotope')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        neutron = initializeNewNeutron()
        neutron._energy = 1.0
        
        try:
            mod.sampleIsotope(neutron)
        except:
            self.fail('Unable to sample isotope')


    # Test Material collideNeutron function
    def testCollideNeutron(self):
        py_printf('UNITTEST', 'Testing Material collideNeutron')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        neutron = initializeNewNeutron()
        neutron._energy = 1.0
        
        try:
            mod.collideNeutron(neutron)
        except:
            self.fail('Unable to collide neutron')


    # Test Material clone function
    def testClone(self):
        py_printf('UNITTEST', 'Testing Material clone')
        mod = Material('mod')
        mod.setDensity(5., 'g/cc')
        mod.addIsotope(self.h1, 1.0)
        mod_clone = mod.clone()
        self.assertEqual(mod_clone.getMaterialNumberDensity(),mod.getMaterialNumberDensity())


class TestRegion(unittest.TestCase):

    def setUp(self):
    
        # Set main simulation params
        self.num_batches = 5
        self.num_neutrons_per_batch = 10000
        self.num_threads = 4
    
        # Define isotopes
        self.h1 = Isotope('H-1')
        self.o16 = Isotope('O-16')
        self.u235 = Isotope('U-235')
        self.u238 = Isotope('U-238')
    
        # Define materials
        self.mix = Material('fuel moderator mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.o16, 1.0)
        self.mix.addIsotope(self.u238, 0.50)
        self.mix.addIsotope(self.u235, .005)

    # Test Region getRegionName function
    def testGetRegionName(self):
        py_printf('UNITTEST', 'Testing Region getRegionName')
        region_fuel = Region('fuel', FUEL)
        self.assertEqual(region_fuel.getRegionName(), 'fuel')


    # Test Region getVolume function
    def testGetVolume(self):
        py_printf('UNITTEST', 'Testing Region getVolume')
        region_fuel = Region('fuel', FUEL)
        region_fuel.setVolume(1.0)
        self.assertAlmostEqual(region_fuel.getVolume(), 1.0)


    # Test Region getMaterial function
    def testGetMaterial(self):
        py_printf('UNITTEST', 'Testing Region getMaterial')
        region_mix = Region('region_mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertEqual(region_mix.getMaterial(), self.mix)


    # Test Region getMaterial function
    def testGetMaterial(self):
        py_printf('UNITTEST', 'Testing Region getMaterial')
        region_mix = Region('region_mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertEqual(region_mix.getMaterial().getMaterialName(), self.mix.getMaterialName())

    # Test Region containsIsotope function
    def testContainsIsotope(self):
        py_printf('UNITTEST', 'Testing Region containsIsotope')
        region_mix = Region('region_mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertTrue(region_mix.containsIsotope(self.h1))


    # Test Region getRegionType function
    def testGetRegionType(self):
        py_printf('UNITTEST', 'Testing Region getRegionType')
        region_mix = Region('region_mix', INFINITE)
        self.assertEqual(region_mix.getRegionType(), 2)

    # Test Region isModerator function
    def testIsModerator(self):
        py_printf('UNITTEST', 'Testing Region isModerator')
        region_mod = Region('mod', MODERATOR)
        self.assertTrue(region_mod.isModerator())

            
    # Test Region isFuel function
    def testIsFuel(self):
        py_printf('UNITTEST', 'Testing Region isFuel')
        region_fuel = Region('fuel', FUEL)
        self.assertTrue(region_fuel.isFuel())


    # Test Region isInfinite function
    def testIsInfinite(self):
        py_printf('UNITTEST', 'Testing Region isInfinite')
        region_mix = Region('mix', INFINITE)
        self.assertTrue(region_mix.isInfinite())


    # Test Region getFuelRadius function
    def testGetFuelRadius(self):
        py_printf('UNITTEST', 'Testing Region getFuelRadius')
        region_fuel = Region('fuel', FUEL)
        region_fuel.setMaterial(self.mix)
        region_fuel.setFuelRadius(1.1)
        self.assertAlmostEqual(region_fuel.getFuelRadius(), 1.1)


    # Test Region getPitch function
    def testGetPitch(self):
        py_printf('UNITTEST', 'Testing Region getPitch')
        region_fuel = Region('mix', FUEL)
        region_fuel.setMaterial(self.mix)
        region_fuel.setPitch(1.1)
        self.assertAlmostEqual(region_fuel.getPitch(), 1.1)


    # Test Region getBucklingSquared function
    def testGetBucklingSquared(self):
        py_printf('UNITTEST', 'Testing Region getBucklingSquared')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        region_mix.setBucklingSquared(1.1)
        self.assertAlmostEqual(region_mix.getBucklingSquared(), 1.1)


    # Test Region getTotalMacroXS function
    def testGetTotalMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Region getTotalMacroXS by Energy')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getTotalMacroXS(1.0), 0.0)

    
    # Test Region getTotalMacroXS function
    def testGetTotalMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Region getTotalMacroXS by Index')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getTotalMacroXS(100), 0.0)


    # Test Region getElasticMacroXS function
    def testGetElasticMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Region getElasticMacroXS by Energy')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getElasticMacroXS(1.0), 0.0)
    
    
    # Test Region getElasticMacroXS function
    def testGetElasticMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Region getElasticMacroXS by Index')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getElasticMacroXS(100), 0.0)


    # Test Region getTotalMacroXS function
    def testGetAbsorptionMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Region getAbsorptionMacroXS by Energy')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getAbsorptionMacroXS(1.0), 0.0)
    
    
    # Test Region getAbsorptionMacroXS function
    def testGetAbsorptionMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Region getAbsorptionMacroXS by Index')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getAbsorptionMacroXS(100), 0.0)


    # Test Region getCaptureMacroXS function
    def testGetCaptureMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Region getCaptureMacroXS by Energy')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getCaptureMacroXS(1.0), 0.0)
    
    
    # Test Region getCaptureMacroXS function
    def testGetCaptureMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Region getCaptureMacroXS by Index')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getCaptureMacroXS(100), 0.0)


    # Test Region getFissionMacroXS function
    def testGetFissionMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Region getFissionMacroXS by Energy')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getFissionMacroXS(1.0), 0.0)
    
    
    # Test Region getFissionMacroXS function
    def testGetFissionMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Region getFissionMacroXS by Index')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getFissionMacroXS(100), 0.0)


    # Test Region getTransportMacroXS function
    def testGetTransportMacroXSbyEnergy(self):
        py_printf('UNITTEST', 'Testing Region getTransportMacroXS by Energy')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getTransportMacroXS(1.0), 0.0)
    
    
    # Test Region getTransportMacroXS function
    def testGetTransportMacroXSbyIndex(self):
        py_printf('UNITTEST', 'Testing Region getTransportMacroXS by Index')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        self.assertGreater(region_mix.getTransportMacroXS(100), 0.0)


    # Test Region collideNeutron function
    def testCollideNeutron(self):
        py_printf('UNITTEST', 'Testing Region collideNeutron')
        region_mix = Region('mix', INFINITE)
        region_mix.setMaterial(self.mix)
        neutron = initializeNewNeutron()
        neutron._energy = 1.0
        
        try:
            region_mix.collideNeutron(neutron)
        except:
            self.fail('Unable to collide neutron')


    # Test Region contains function
    def testContains(self):
        py_printf('UNITTEST', 'Testing Region contains')
        region_fuel = Region('mix', FUEL)
        region_fuel.setMaterial(self.mix)
        region_fuel.setFuelRadius(1.0)
        region_fuel.setFuelRadius(2.5)
        neutron = initializeNewNeutron()
        neutron._x = 0.5
        neutron._y = 0.5
        self.assertTrue(region_fuel.contains(neutron))


    # Test Region onBoundary function
    def testOnBoundary(self):
        py_printf('UNITTEST', 'Testing Region onBoundary')
        region_fuel= Region('fuel', FUEL)
        region_fuel.setFuelRadius(0.5)
        region_fuel.setPitch(1.5)
        neutron = initializeNewNeutron()
        neutron._x = 0.5
        neutron._y = 0.0
        self.assertTrue(region_fuel.onBoundary(neutron))



class TestGeometry(unittest.TestCase):

    def setUp(self):
    
        # Set main simulation params
        self.num_batches = 10
        self.num_neutrons_per_batch = 10000
        self.num_threads = 5
    
        # Define isotopes
        self.h1 = Isotope('H-1')
        self.o16 = Isotope('O-16')
        self.u235 = Isotope('U-235')
        self.u238 = Isotope('U-238')
        self.b10 = Isotope('B-10')
        self.zr90 = Isotope('Zr-90')
        
        # Define materials
        self.mix = Material('mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.o16, 1.0)
        self.mix.addIsotope(self.u238, 0.50)
        self.mix.addIsotope(self.u235, .025)
        self.mix.addIsotope(self.zr90, .16)
        self.mix.addIsotope(self.b10, .0000001)
    
        # Define region
        self.region_mix = Region('mix', INFINITE)
        self.region_mix.setMaterial(self.mix)
    
    
    # Test Geometry getNumNeutronsPerBatch function
    def testGetNumNeutronsPerBatch(self):
        py_printf('UNITTEST', 'Testing Geometry getNumNeutronsPerBatch')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        self.assertEqual(self.num_neutrons_per_batch, geometry.getNumNeutronsPerBatch())
    
    
    # Test Geometry getTotalNumNeutrons function
    def testGetTotalNumNeutrons(self):
        py_printf('UNITTEST', 'Testing Geometry getTotalNumNeutrons')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        geometry.setNumBatches(self.num_batches)
        self.assertEqual(100000, geometry.getTotalNumNeutrons())


    # Test Geometry getNumBatches function
    def testGetNumBatches(self):
        py_printf('UNITTEST', 'Testing Geometry getNumBatches')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.setNumBatches(self.num_batches)
        self.assertEqual(10, geometry.getNumBatches())


    # Test Geometry getNumThreads function
    def testGetNumThreads(self):
        py_printf('UNITTEST', 'Testing Geometry getNumThreads')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.setNumThreads(self.num_threads)
        self.assertEqual(5, geometry.getNumThreads())


    # Test Geometry getSpatialType function
    def testGetSpatialType(self):
        py_printf('UNITTEST', 'Testing Geometry getSpatialType')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        self.assertEqual(0, geometry.getSpatialType())


    # Test Geometry getBucklingSquared function
    def testGetBucklingSquared(self):
        py_printf('UNITTEST', 'Testing Geometry getBucklingSquared')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.setBucklingSquared(1.1)
        self.assertAlmostEqual(1.1, geometry.getBucklingSquared())


    # Test Geometry getVolume function
    def testGetVolume(self):
        py_printf('UNITTEST', 'Testing Geometry getVolume')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.addRegion(self.region_mix)
        self.assertAlmostEqual(self.region_mix.getVolume(), geometry.getVolume())


    # Test Geometry setDancoff function
    def testSetDancoff(self):
        py_printf('UNITTEST', 'Testing Geometry setDancoff')
        geometry = Geometry(HOMOGENEOUS_EQUIVALENCE)
        region_fuel = Region('fuel', FUEL)
        region_fuel.setMaterial(self.mix)
        region_mod = Region('mod', MODERATOR)
        region_mod.setMaterial(self.mix)
        geometry.addRegion(region_mod)
        geometry.addRegion(region_fuel)

        try:
            geometry.setDancoffFactor(0.7)
        except:
            self.fail('Could not set dancoff number')


    # Test Geometry addRegion function
    def testAddRegion(self):
        py_printf('UNITTEST', 'Testing Geometry addRegion')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        
        try:
            geometry.addRegion(self.region_mix)
        except:
            self.fail('Could not add Region')


    # Test Geometry runMonteCarloSimulation function
    def testRunMonteCarloSimulation(self):
        py_printf('UNITTEST', 'Testing Geometry runMonteCarloSimulation')
        geometry = Geometry(INFINITE_HOMOGENEOUS)
        geometry.addRegion(self.region_mix)
        geometry.setNumBatches(self.num_batches)
        geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        geometry.setNumThreads(self.num_threads)
    
        try:
            geometry.runMonteCarloSimulation()
        except:
            self.fail('Could not run Monte Carlo')


#class TestProcess(unittest.TestCase):
#
#    def setUp(self):
#    
#    
#    
#    def testGetTallyCenters(self):
#    
#    def testGetTallyEdges(self):
#    
#    def testGetTallyBatchMu(self):
#    
#    def testGetTallyBatchVariances(self):
#    
#    def testGetTallyBatchStdDev(self):
#    
#    def testGetTallyBatchRelError(self):
#    
#    def testGetBatchTallyStatistics(self):
#    
#    def testPrintTallies(self):
#    
#    def testComputeMeanNumCollisions(self):
#    
#    def testComputeMeanNeutronLifetime(self):
#
#
#    def testRIEffSetName(self):
#
#    def testRIEffComputeRIs(self):
#
#    def testRIEffGetNumIntegrals(self):
#
#    def testRIEffGetIntegrals(self):
#
#    def testRIEffGetVariances(self):
#
#    def testRIEffGetStandardDeviations(self):
#
#    def testRIEffGetRelativeError(self):
#
#    def testRIEffGetEnergyBandsCenters(self):
#
#    def testRIEffGetEnergyBands(self):
#
#    def testRIEffGetRITally(self):
#
#    def testRIEffprintRI(self):
#
#    def testRIEffOutputRItoFile(self):
#
#
#    def testRITrueSetName(self):
#
#    def testRITrueComputeRIs(self):
#
#    def testRITrueGetNumIntegrals(self):
#    
#    def testRITrueGetIntegrals(self):
#    
#    def testRITrueGetVariances(self):
#    
#    def testRITrueGetStandardDeviations(self):
#    
#    def testRITrueGetRelativeError(self):
#    
#    def testRITrueGetEnergyBandsCenters(self):
#    
#    def testRITrueGetEnergyBands(self):
#    
#    def testRITrueGetRITally(self):
#    
#    def testRITrueprintRI(self):
#    
#    def testRITrueOutputRItoFile(self):
#
#    def testPrintRIs(self):
#
#    
#    def testGroupXSSetName(self):
#    
#    def testGroupXSComputeRIs(self):
#    
#    def testGroupXSGetNumXS(self):
#    
#    def testGroupXSGetXS(self):
#    
#    def testGroupXSGetVariances(self):
#    
#    def testGroupXSGetStandardDeviations(self):
#    
#    def testGroupXSGetRelativeError(self):
#    
#    def testGroupXSGetEnergyBandsCenters(self):
#    
#    def testGroupXSGetEnergyBands(self):
#    
#    def testGroupXSGetRITally(self):
#    
#    def testGroupXSprintXS(self):
#    
#    def testGroupXSOutputXStoFile(self):





#class TestPlotting(unittest.TestCase)
#class TestSLBW(unittest.TestCase)
#class TestSurface(unittest.TestCase)







