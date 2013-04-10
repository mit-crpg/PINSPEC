import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *
import unittest
import os


class TestInfinite(unittest.TestCase):
    
    def setUp(self):

        # Set main simulation params
        self.num_batches = 10
        self.num_neutrons_per_batch = 10000
        self.num_threads = 5
        self.nu = 2.45

        # Define isotopes
        self.h1 = Isotope('H-1')
        self.o16 = Isotope('O-16')
        self.u235 = Isotope('U-235')
        self.u238 = Isotope('U-238')
        self.zr90 = Isotope('Zr-90')
        self.b10 = Isotope('B-10')
    
    
    # check Keff for infinite geometry
    def testInfiniteKeff(self):

        py_printf('UNITTEST', 'Testing Infinite keff')
        
        # Define materials
        self.mix = Material('fuel moderator mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.b10, .0000001)
        self.mix.addIsotope(self.o16, 1.0)
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.u238, 0.5)
        self.mix.addIsotope(self.u235, 0.025)
        self.mix.addIsotope(self.zr90, 0.16)
        
        # Define regions
        self.region_mix = Region('infinite medium', INFINITE)
        self.region_mix.setMaterial(self.mix)
        
        # Define geometry
        self.geometry = Geometry(INFINITE_HOMOGENEOUS)
        self.geometry.addRegion(self.region_mix)
        
        # Create a tally for the flux
        self.fission_tally = createTally(self.region_mix, FISSION_RATE)
        self.abs_tally = createTally(self.region_mix, ABSORPTION_RATE)
        self.fission_tally.generateBinEdges(1E-2, 1E7, 1, EQUAL)
        self.abs_tally.generateBinEdges(1E-2, 1E7, 1, EQUAL)
        TallyBank.registerTally(self.fission_tally)
        TallyBank.registerTally(self.abs_tally)
        
        # Run Monte Carlo simulation
        self.geometry.runMonteCarloSimulation()
        
        self.fission_tally.computeBatchStatistics()
        self.abs_tally.computeBatchStatistics()
    
        k_eff_tally = self.fission_tally / self.abs_tally
        k_eff_array = k_eff_tally.retrieveTallyMu(1)
        k_eff = self.nu * k_eff_array[0]

        self.assertGreater(.05, abs(k_eff - 1.2721)/1.2721)

    
    # check 6-10 eV resonance integral for infinite dilute case
    def testInfiniteDiluteRI6_10(self):

        py_printf('UNITTEST', 'Testing Infinite RI 6-10 eV')
        
        # Define materials
        self.mix = Material('fuel moderator mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.u238, 1E-6)
        
        # Define regions
        self.region_mix = Region('infinite medium', INFINITE)
        self.region_mix.setMaterial(self.mix)
        
        # Define geometry
        self.geometry = Geometry(INFINITE_HOMOGENEOUS)
        self.geometry.addRegion(self.region_mix)
        self.geometry.setNumBatches(self.num_batches)
        self.geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        self.geometry.setNumThreads(self.num_threads)
        
        # Create a tally for the flux
        u238_abs_rate = createTally(self.u238, ABSORPTION_RATE)
        flux = createTally(self.region_mix, FLUX)
        
        abs_rate_bin_edges = numpy.array([1E-5, 1., 6., 10., 25., 50., 100., 1000.])
        u238_abs_rate.setBinEdges(abs_rate_bin_edges)
        flux.setBinEdges(abs_rate_bin_edges)
        
        TallyBank.registerTally(u238_abs_rate, self.region_mix)
        TallyBank.registerTally(flux)

        # Run Monte Carlo simulation
        self.geometry.runMonteCarloSimulation()
        
        RI_eff = process.RIEff(flux, u238_abs_rate, 'RI infinite')
        RI_eff.printRI()
        
        RIs = RI_eff._RI.retrieveTallyMu(8)
        
        self.assertGreater(.05, abs(127.7-RIs[2])/127.7)


    # check 10-25 eV resonance integral for infinite dilute case
    def testInfiniteDiluteRI10_25(self):
        
        py_printf('UNITTEST', 'Testing Infinite RI 10-25 eV')
        
        # Define materials
        self.mix = Material('fuel moderator mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.u238, 1E-6)
        
        # Define regions
        self.region_mix = Region('infinite medium', INFINITE)
        self.region_mix.setMaterial(self.mix)
        
        # Define geometry
        self.geometry = Geometry(INFINITE_HOMOGENEOUS)
        self.geometry.addRegion(self.region_mix)
        self.geometry.setNumBatches(self.num_batches)
        self.geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        self.geometry.setNumThreads(self.num_threads)
        
        # Create a tally for the flux
        u238_abs_rate = createTally(self.u238, ABSORPTION_RATE)
        flux = createTally(self.region_mix, FLUX)
        
        abs_rate_bin_edges = numpy.array([1E-5, 1., 6., 10., 25., 50., 100., 1000.])
        u238_abs_rate.setBinEdges(abs_rate_bin_edges)
        flux.setBinEdges(abs_rate_bin_edges)
        
        TallyBank.registerTally(u238_abs_rate, self.region_mix)
        TallyBank.registerTally(flux)
        
        # Run Monte Carlo simulation
        self.geometry.runMonteCarloSimulation()
        
        RI_eff = process.RIEff(flux, u238_abs_rate, 'RI infinite')
        RI_eff.printRI()
        
        RIs = RI_eff._RI.retrieveTallyMu(8)
        
        self.assertGreater(.05, abs(66.1-RIs[3])/66.1)

    
    # check 25-50 eV resonance integral for infinite dilute case
    def testInfiniteDiluteRI25_50(self):
        
        py_printf('UNITTEST', 'Testing Infinite RI 25-50 eV')
        
        # Define materials
        self.mix = Material('fuel moderator mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.u238, 1E-6)
        
        # Define regions
        self.region_mix = Region('infinite medium', INFINITE)
        self.region_mix.setMaterial(self.mix)
        
        # Define geometry
        self.geometry = Geometry(INFINITE_HOMOGENEOUS)
        self.geometry.addRegion(self.region_mix)
        self.geometry.setNumBatches(self.num_batches)
        self.geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        self.geometry.setNumThreads(self.num_threads)
        
        # Create a tally for the flux
        u238_abs_rate = createTally(self.u238, ABSORPTION_RATE)
        flux = createTally(self.region_mix, FLUX)
        
        abs_rate_bin_edges = numpy.array([1E-5, 1., 6., 10., 25., 50., 100., 1000.])
        u238_abs_rate.setBinEdges(abs_rate_bin_edges)
        flux.setBinEdges(abs_rate_bin_edges)
        
        TallyBank.registerTally(u238_abs_rate, self.region_mix)
        TallyBank.registerTally(flux)
        
        # Run Monte Carlo simulation
        self.geometry.runMonteCarloSimulation()
        
        RI_eff = process.RIEff(flux, u238_abs_rate, 'RI infinite')
        RI_eff.printRI()
        
        RIs = RI_eff._RI.retrieveTallyMu(8)
        
        self.assertGreater(.05, abs(41.8-RIs[4])/41.8)



class TestEquivalence(unittest.TestCase):
    
    def setUp(self):
        
        # Set main simulation params
        self.num_batches = 10
        self.num_neutrons_per_batch = 10000
        self.num_threads = 5
        self.nu = 2.45
        self.radius_fuel = 0.4096
        self.pitch = 1.26
        self.dancoff = 0.277
        
        # Define isotopes
        self.h1 = Isotope('H-1')
        self.o16 = Isotope('O-16')
        self.u235 = Isotope('U-235')
        self.u238 = Isotope('U-238')
        self.zr90 = Isotope('Zr-90')
        self.b10 = Isotope('B-10')
    
    
    # check Keff for equivalence geometry
    def testEquivalenceKeff(self):
        
        py_printf('UNITTEST', 'Testing Equivalence keff')

        # Define materials
        self.fuel = Material('fuel')
        self.mod = Material('moderator')
        self.fuel.setDensity(9., 'g/cc')
        self.mod.setDensity(1., 'g/cc')
        self.mod.addIsotope(self.b10, .0000001)
        self.mod.addIsotope(self.o16, 1.0)
        self.mod.addIsotope(self.h1, 2.0)
        self.fuel.addIsotope(self.o16, 1.05)
        self.fuel.addIsotope(self.u238, 0.5)
        self.fuel.addIsotope(self.u235, 0.025)
        self.fuel.addIsotope(self.zr90, 0.16)
        
        # Define regions
        self.region_mod = Region('moderator', MODERATOR)
        self.region_fuel = Region('fuel', FUEL)
        self.region_mod.setMaterial(self.mod)
        self.region_fuel.setMaterial(self.fuel)
        self.region_fuel.setFuelRadius(self.radius_fuel)
        self.region_fuel.setPitch(self.pitch)
        self.region_mod.setFuelRadius(self.radius_fuel)
        self.region_mod.setPitch(self.pitch)
        
        # Define geometry
        self.geometry = Geometry(HOMOGENEOUS_EQUIVALENCE)
        self.geometry.addRegion(self.region_mod)
        self.geometry.addRegion(self.region_fuel)
        self.geometry.setNumBatches(self.num_batches)
        self.geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        self.geometry.setNumThreads(self.num_threads)
        self.geometry.setDancoffFactor(self.dancoff)
        
        # Create a tally for the flux
        self.fission_tally = createTally(self.geometry, FISSION_RATE)
        self.abs_tally = createTally(self.geometry, ABSORPTION_RATE)
        self.fission_tally.generateBinEdges(1E-2, 1E7, 1, EQUAL)
        self.abs_tally.generateBinEdges(1E-2, 1E7, 1, EQUAL)
        TallyBank.registerTally(self.fission_tally)
        TallyBank.registerTally(self.abs_tally)
        
        # Run Monte Carlo simulation
        self.geometry.runMonteCarloSimulation()
        
        self.fission_tally.computeBatchStatistics()
        self.abs_tally.computeBatchStatistics()
        
        k_eff_tally = self.fission_tally / self.abs_tally
        k_eff_array = k_eff_tally.retrieveTallyMu(1)
        k_eff = self.nu * k_eff_array[0]
        
        self.assertGreater(.05, abs(k_eff - 1.5623)/1.5623)


