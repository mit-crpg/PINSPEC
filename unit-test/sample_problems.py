import numpy
from pinspec import *
import pinspec.SLBW as SLBW
import pinspec.plotter as plotter
import pinspec.process as process
from pinspec.log import *
import math
import unittest


# run 2012 22.211 hw2
class TestHW2(unittest.TestCase):
    
    def setUp(self):
    
        # Set main simulation params
        self.num_neutrons = 100000
        setOutputDirectory('HW2')
        
        # Initialize isotopes
        h1 = Isotope('H-1')
        c12 = Isotope('C-12')
        
        # Create an artifical capture xs for hydrogen
        norm_const = 0.025
        h1_capture_energies = numpy.logspace(-5., 7.5, 500);
        h1_capture_xs = numpy.sqrt(norm_const/h1_capture_energies) * 7.
        h1.setCaptureXS(h1_capture_energies, h1_capture_xs)
        
        
        ###########################################################################
        ###########################   Problems 3-6   ##############################
        ###########################################################################
        
        h1_material = Material('H-1')
        h1_material.setDensity(0.07778, 'g/cc')
        h1_material.addIsotope(h1, 2.0)
        
        flux = createTally(h1_material, FLUX)
        flux.generateBinEdges(1E-2, 1E7, 1000, LOGARITHMIC)
        
        self.coll_rate_1eV = createTally(h1_material, COLLISION_RATE)
        self.coll_rate_1eV.generateBinEdges(1E-1, 2E6, 1, EQUAL)
        
        self.coll_rate = createTally(h1_material, COLLISION_RATE)
        self.coll_rate.generateBinEdges(1E-7, 2E6, 1, EQUAL)
        
        self.times = createTally(h1_material, INTERCOLLISION_TIME)
        self.times.generateBinEdges(1E-7, 2E6, 1, EQUAL)
        
        fissioner = Fissioner()
        neutron = initializeNewNeutron()
        
        for i in range(self.num_neutrons):
            
            # Sample a fission energy from the Watt fission spectrum
            neutron._energy = fissioner.emitNeutroneV()
            neutron._alive = True
            reached_one_ev = False
            
            # Simulate neutron until it is absorbed in H-1
            while(neutron._alive):
                
                h1_material.collideNeutron(neutron)
                flux.tally(neutron)
                self.times.tally(neutron)
                self.coll_rate.tally(neutron)
                
                if neutron._energy < 1.0:
                    reached_one_ev = True
                
                if not reached_one_ev:
                    self.coll_rate_1eV.tally(neutron)


    
    def testMeanCollToOne(self):
        
        py_printf('UNITTEST', 'Testing HW2 mean collisions to 1 eV')
        
        num_collisions = process.computeMeanNumCollisions(self.coll_rate_1eV, \
                                                      self.num_neutrons)
                        
        self.assertGreater(.05, abs(13.6 - num_collisions)/13.6)

    def testMeanColl(self):
        
        py_printf('UNITTEST', 'Testing HW2 mean collisions')
    
        num_collisions = process.computeMeanNumCollisions(self.coll_rate, \
                                                      self.num_neutrons)
        
        self.assertGreater(.05, abs(21.0 - num_collisions)/21.0)
    
    
    def testMeanLifetime(self):
        
        py_printf('UNITTEST', 'Testing HW2 mean lifetime')
        
        mean_lifetime = process.computeMeanNeutronLifetime(self.times, self.num_neutrons)
        
        self.assertGreater(.10, abs(1.40E-5 - mean_lifetime) / 1.40E-5)

