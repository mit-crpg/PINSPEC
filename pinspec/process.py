import numpy as np
from pinspec import *
import matplotlib.pyplot as plt
from log import *


class RI(object):

	def __init__(self, flux, abs_rate):
	
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
			py_printf('NORMAL', 'RI [ %9.4f - %9.4f eV] =  %9.4f', self.bin_edges[i], self.bin_edges[i+1], self.RIs[i])

	def getBinCenters(self):
		return self.bin_centers			
			
	def getBinEdges(self):
		return self.bin_edges
	
	def getRIs(self):
		return self.RIs
	
	

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
				py_printf('NORMAL', 'Group XS [ %9.4f - %9.4f eV] =  %9.3E cm^-1', self.bin_edges[i], self.bin_edges[i+1], self.groupXS[i,j])

