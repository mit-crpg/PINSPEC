import numpy as np
import matplotlib.pyplot as plt
from pinspec import *
from log import *


# NEEDS LOTS OF ERROR CHECKING!!!!

class ReactionRate(object):

    def __init__(self, reaction_rate):
    
        self.reaction_rate = reaction_rate
        self.num_rates = reaction_rate.getNumBins()
        self.reaction_rate_mu = reaction_rate.retrieveTallyMu(self.num_rates)
        self.bin_edges = reaction_rate.retrieveTallyEdges(self.num_rates+1)
        self.bin_centers = reaction_rate.retrieveTallyCenters(self.num_rates)
                
    def printReactionRates(self):

        for i in range(0,self.num_rates):
            py_printf('RESULT', 'Reaction Rate [ %7.2f - %7.2f eV  ] =  %5.2f', \
                            self.bin_edges[i], self.bin_edges[i+1], \
                            self.reaction_rate_mu[i])

    def getNumRates(self):
        return self.num_rates

    def getBinCenters(self):
        return self.bin_centers            
            
    def getBinEdges(self):
        return self.bin_edges
    
    def getReactionRates(self):
        return self.reaction_rate_mu

    def __div__(self, other_rate):

        return ReactionRateRatio(self, other_rate)


class ReactionRateRatio(object):

    def __init__(self, rate1, rate2):

        self.name = ''
        self.rate1 = rate1
        self.rate2 = rate2
        self.num_rates1 = rate1.getNumRates()
        self.num_rates2 = rate2.getNumRates()
        self.bin_edges1 = rate1.getBinEdges()
        self.bin_edges2 = rate2.getBinEdges()
        self.bin_centers1 = rate1.getBinCenters()
        self.bin_centers2 = rate2.getBinCenters()

        if self.num_rates2 is not 1:
            if self.num_rates1 is not self.num_rates2:
                py_printf('ERROR', 'Unable to divide two reaction rates' \
                                        ' with different numbers of rates')
            elif not (self.bin_edges1==self.bin_edges2).all():
                py_printf('ERROR', 'Unable to divide two reaction rates' \
                                                    ' with different bin edges')
        else:
            self.ratios = rate1.getReactionRates() / rate2.getReactionRates()

                
    def printRatios(self):

        for i in range(0,self.num_rates1):
            py_printf('RESULT', 'Reaction Rate Ratio [ %7.2f - %7.2f eV  ] =  '\
                        ' %5.4f', self.bin_edges1[i], self.bin_edges1[i+1], \
                                                                self.ratios[i])

    def getNumRatios(self):
        return self.num_rates1

    def getBinCenters(self):
        return self.bin_centers1        
            
    def getBinEdges(self):
        return self.bin_edges1
    
    def getRatios(self):
        return self.ratios

    def setName(self, name):
        self.name = str(name)

    def getName(self):
        return self.name


def printReactionRateRatios(ratios):

        num_bins = ratios[0].getNumRatios()
        bin_edges = ratios[0].getBinEdges()
        bin_centers = ratios[0].getBinCenters()

        title = 'Printing reaction rate ratios...'.ljust(50)
        for ratio in ratios:
            title += ratio.getName().center(12)

        py_printf('RESULT', title)

        for i in range(num_bins):
            string = 'Reaction Rate Ratio [ %7.2f - %7.2f eV  ] =  ' % \
                                                  (bin_edges[i], bin_edges[i+1])
            for ratio in ratios:
                string += ('%4.4f' % ratio.getRatios()[i]).center(12)

            py_printf('RESULT', string)


def printResonanceIntegrals(integrals):

        num_bins = integrals[0].getNumIntegrals()
        bin_edges = integrals[0].getBinEdges()
        bin_centers = integrals[0].getBinCenters()

        title = 'Printing resonance integrals...'.ljust(48)
        for integral in integrals:
            title += integral.getName().center(12)

        py_printf('RESULT', title)

        for i in range(num_bins):
            string = 'Resonance Integral [ %7.2f - %7.2f eV  ] =  ' % \
                                                  (bin_edges[i], bin_edges[i+1])
            for integral in integrals:
                string += ('%4.4f' % integral.getRIs()[i]).center(12)

            py_printf('RESULT', string)
    

class RI(object):

    def __init__(self, flux, abs_rate):
    
        self.name = ''
        self.num_RIs = abs_rate.getNumBins()
        self.flux = flux
        self.abs_rate = abs_rate
        self.flux_mu = flux.retrieveTallyMu(self.num_RIs)
        self.abs_rate_mu = abs_rate.retrieveTallyMu(self.num_RIs)
        self.bin_edges = abs_rate.retrieveTallyEdges(self.num_RIs+1)
        self.bin_centers = abs_rate.retrieveTallyCenters(self.num_RIs)

        self.RIs = np.zeros(self.num_RIs)
    
        self.computeRI()
        
        
    def computeRI(self):
        
        for i in range(0,self.num_RIs):
            self.RIs[i] = self.abs_rate_mu[i] / self.flux_mu[i] * np.log(self.bin_edges[i+1]/self.bin_edges[i])
                
                
    def printRI(self):

        for i in range(0,self.num_RIs):
            py_printf('RESULT', 'Resonance Integral [ %7.2f - %7.2f eV  ] =  ' \
                                            '%5.2f', self.bin_edges[i], 
                                            self.bin_edges[i+1], self.RIs[i])
    def getNumIntegrals(self):
        return self.num_RIs

    def getBinCenters(self):
        return self.bin_centers            
            
    def getBinEdges(self):
        return self.bin_edges
    
    def getRIs(self):
        return self.RIs

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name
    
    

class groupXS(object):
    
    def __init__(self, flux, rxns):
        
        self.flux = flux    
        
        if type(rxns) is not list:            
            rxns = [rxns]
        
        self.num_bins = rxns[0].getNumBins()
        self.bin_edges = rxns[0].retrieveTallyEdges(self.num_bins+1)
        self.rxns = rxns
        self.num_rxns = len(rxns)
        
        self.groupXS = np.zeros(shape=(self.num_bins, self.num_rxns))
        self.computeGroupXS()
        
    def computeGroupXS(self):
        
        # Get the bin center energies
        n_bins = self.flux.getNumBins()
        self.flux.computeBatchStatistics()
        flux_bin_centers = self.flux.retrieveTallyCenters(n_bins)

        flux_mu = self.flux.retrieveTallyMu(n_bins)
        for j in range(0,self.num_rxns):
            rxn_mu = self.rxns[j].retrieveTallyMu(self.num_bins) 
            for i in range(0,self.num_bins):
                self.groupXS[i,j] = rxn_mu[i] / flux_mu[i]
                
                
                
    def printGroupXS(self):
        
        for j in range(0,self.num_rxns):
            rxn_mu = self.rxns[j].retrieveTallyMu(self.num_bins) 
            for i in range(0,self.num_bins):
                py_printf('RESULT', 'Group XS [ %7.2f - %7.2f eV  ] =  %5.3E cm^-1', self.bin_edges[i], self.bin_edges[i+1], self.groupXS[i,j])



def computeMeanNumCollisions(coll_rate, num_neutrons):

    coll_rate.computeScaledBatchStatistics(num_neutrons)

    num_bins = coll_rate.getNumBins()
    coll_rate_mu = coll_rate.retrieveTallyMu(num_bins)
    mean_rate = 0.0

    for i in range(num_bins):
        mean_rate += coll_rate_mu[i]

    return mean_rate


def computeMeanNeutronLifetime(coll_times, num_neutrons):

    coll_times.computeScaledBatchStatistics(num_neutrons)

    num_bins = coll_times.getNumBins()
    coll_times_mu = coll_times.retrieveTallyMu(num_bins)
    mean_time = 0.0

    for i in range(num_bins):
        mean_time += coll_times_mu[i]

    return mean_time

