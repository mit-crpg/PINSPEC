from pinspec import *
import unittest
import numpy

class TestInfiniteMedium(unittest.TestCase):

    def setUp(self):

        # Set main simulation params
        self.num_batches = 5
        self.num_neutrons_per_batch = 10000
        self.num_threads = 4
        log_setlevel(INFO)

        setXSLibDirectory('../xs-lib/')   # This is also a default, but set it as example

        # Define isotopes
        self.h1 = Isotope('H-1')
        self.o16 = Isotope('O-16')
        self.u235 = Isotope('U-235')
        self.u238 = Isotope('U-238')
        
        # Define materials
        self.mix = Material()
        self.mix.setMaterialName('fuel moderator mix')
        self.mix.setDensity(5., 'g/cc')
        self.mix.addIsotope(self.h1, 1.0)
        self.mix.addIsotope(self.o16, 1.0)
        self.mix.addIsotope(self.u238, 0.50)
        self.mix.addIsotope(self.u235, .005)

        # Define regions
        self.region_mix = Region('infinite medium', INFINITE)
        self.region_mix.setMaterial(self.mix)


    def testTallyNumLogarithmicBins(self):
        flux = Tally('flux test', GEOMETRY, FLUX)
        flux.generateBinEdges(1E-2, 1E7, 2000, LOGARITHMIC)
        self.assertEqual(2000, flux.getNumBins())


    def testTallyNumEqualBins(self):
        flux = Tally('flux test', GEOMETRY, FLUX)
        flux.generateBinEdges(1E-2, 1E7, 2000, EQUAL)
        self.assertEqual(2000, flux.getNumBins())


    def testTallyNumUserDefinedBins(self):
        abs_rate = Tally('flux test', GEOMETRY, ABSORPTION_RATE)
        bin_edges = numpy.array([1.0, 2.25, 3.0, 100.5])
        abs_rate.setBinEdges(bin_edges)
        self.assertEqual(3, abs_rate.getNumBins())


    def testUnAddedRegion(self):
        geometry = Geometry()
        geometry.setSpatialType(INFINITE_HOMOGENEOUS)
        geometry.setNumBatches(self.num_batches)
        geometry.setNeutronsPerBatch(self.num_neutrons_per_batch)
        geometry.setNumThreads(self.num_threads)
        self.assertRaises(Exception, geometry.runMonteCarloSimulation)


    def testSetBatchNum(self):
        geometry = Geometry()
        geometry.setSpatialType(INFINITE_HOMOGENEOUS)
        geometry.setNumBatches(self.num_batches)
        self.assertEqual(self.num_batches, geometry.getNumBatches())


if __name__ == '__main__':
    unittest.main()

