import numpy as np
import matplotlib.pyplot as plt
from pinspec import *
from log import *
import scipy.integrate as integrate


# NEEDS LOTS OF ERROR CHECKING AND PYTHON DOCSTRING COMMENTING!!!!


###############################################################################
########################  Tally Data Retrieval Methods  #######################
###############################################################################

def getTallyCenters(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally centers from input of type' \
                                                            + str(type(tally)))
    else:
        num_bins = tally.getNumBins()
        centers = tally.retrieveTallyCenters(num_bins)
        return centers


def getTallyEdges(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally edges from input of type' \
                                                            + str(type(tally)))
    else:
        num_bins = tally.getNumBins()
        edges = tally.retrieveTallyEdges(num_bins+1)
        return edges


def getTallyBatchMu(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally mu from input of type' \
                                                        + str(type(tally)))
    if not tally.hasComputedBatchStatistics():
        py_printf('WARNING', 'Unable to get tally mu for tally %s since it' \
                       + ' has not yet computed batch statistics', \
                       + tally.getTallyName())
    else:
        num_bins = tally.getNumBins()
        mu = tally.retrieveTallyMu(num_bins)
        return mu
    

def getTallyBatchVariances(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally variances from input of ' \
                                                  + ' type' + str(type(tally)))
    if not tally.hasComputedBatchStatistics():
        py_printf('WARNING', 'Unable to get tally variances for tally %s'
                        ' since it has not yet computed batch statistics', \
                       + tally.getTallyName())
    else:
        num_bins = tally.getNumBins()
        variances = tally.retrieveTallyVariance(num_bins)
        return variances


def getTallyBatchStdDev(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally std. dev. from input of ' \
                                                  + ' type' + str(type(tally)))
    if not tally.hasComputedBatchStatistics():
        py_printf('WARNING', 'Unable to get tally std. dev. for tally %s'
                        ' since it has not yet computed batch statistics', \
                       + tally.getTallyName())
    else:
        num_bins = tally.getNumBins()
        std_dev = tally.retrieveTallyStdDev(num_bins)
        return std_dev


def getTallyBatchRelErr(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally rel. err. from input of ' \
                                                 + ' type' + str(type(tally)))
    if not tally.hasComputedBatchStatistics():
        py_printf('WARNING', 'Unable to get tally rel. err. for tally %s'
                        ' since it has not yet computed batch statistics', \
                       + tally.getTallyName())
    else:
        num_bins = tally.getNumBins()
        rel_err = tally.retrieveTallyRelErr(num_bins)
        return rel_err


def getBatchTallyStatistics(tally):

    if not isinstance(tally, Tally):
        py_printf('WARNING', 'Unable to get tally statistics from input of ' \
                                                + ' type' + str(type(tally)))
    if not tally.hasComputedBatchStatistics():
        py_printf('WARNING', 'Unable to get tally statistics for tally %s'
                        ' since it has not yet computed batch statistics', \
                       + tally.getTallyName())
    else:
        edges = getTallyEdges(tally)
        centers = getTallyCenters(tally)
        mu = getTallyBatchMu(tally)
        variances = getTallyBatchVariances(tally)
        std_dev = getTallyBatchStdDev(tally)
        rel_err = getTallyBatchRelErr(tally)

        # Create a 2D array of all arrays        
        statistics = np.array([edges, centers, mu, variances, std_dev, rel_err])
        return statistics


#NOTE: This is used to print collections (lists) of tally objects. Since
# RIEff objects are simply Python wrappers for an underlying DERIVED type
# tally, this method can also be used to print lists of RIEff objects. It also
# works for RITrue for now, though RITrue objects are not stored as Tallies
def printTallies(tallies, header='', types='Tallies'):

    if type(tallies) is list:

        # Check that all elements are Tally ojbects
        for tally in tallies:

            if isinstance(tally, RIEff):
                for i in range(len(tallies)):
                    tallies[i] = tallies[i].getRITally()

            elif not isinstance(tally, Tally) and not isinstance(tally, RITrue):
                py_printf('ERROR', 'Unable to print %s since input contains' + \
                        'contains elements of type %s', types, str(type(tally)))

        if isinstance(tallies[0], Tally):

            # Check that all tallies have the same number of bins
            num_bins = tallies[0].getNumBins()
            for tally in tallies:
                if tally.getNumBins() is not num_bins:
                    py_printf('ERROR', 'Unable to print %s since tally %s' \
                                ' has %d bins while tally % has %d bins', types, \
                                tallies[0].getTallyName(), tallies[0].getNumBins(),
                                tally.getTallyName(), tally.getNumBins())

            # Check that all tallies have the same bin edges
            edges = getTallyEdges(tallies[0])
            for tally in tallies:
                if not (edges==getTallyEdges(tally)).all():
                    py_printf('ERROR', 'Unable to print %s since tally %s' \
                                ' has different bin edges than tally %s', types, \
                                tallies[0].getTallyName(), tally.getTallyName())

            # If we passed all checks, then print the tally batch averages        
            # Get the batch mu for each tally
            mu = []
            for tally in tallies:
                mu.append(getTallyBatchMu(tally))


        elif isinstance(tallies[0], RITrue):

            # Check that all RIs have the same number of integrals
            num_bins = tallies[0].getNumIntegrals()
            for tally in tallies:
                if tally.getNumIntegrals() is not num_bins:
                    py_printf('ERROR', 'Unable to print %s since RI %s has %d' \
                            ' integrals while RI % has %d integrals', types, \
                            tallies[0].getName(), tallies[0].getNumIntegrals(),
                            tally.getName(), tally.getNumIntegrals())

            # Check that all RIs have the same bin edges
            edges = tallies[0].getEnergyBands()
            for tally in tallies:
                if not (edges==tally.getEnergyBands()).all():
                    py_printf('ERROR', 'Unable to print %s since RI %s' \
                                ' has different bin edges than RI %s', types, \
                                tallies[0].getName(), tally.getName())

            # If we passed all checks, then print the RIs
            mu = []
            for tally in tallies:
                mu.append(tally.getIntegrals())


        py_printf('HEADER', 'Batch Statistics for ' + header + ' ' + types)
        py_printf('SEPARATOR', '')

        title = ' '.center(7) + 'Energy Band' + ' '.center(9)
        for tally in tallies:
            title += '    Mu    '

        py_printf('RESULT', title)
        py_printf('SEPARATOR', '')

        for i in range(num_bins):
            entry = '[ %7.2f - %7.2f eV ]:' % (edges[i], edges[i+1]) 

            for j in range(len(tallies)):
                if (mu[j][i] <= 1E-2):
                    entry += '  %8.2E' % mu[j][i]
                elif (mu[j][i] >= 10. and mu[j][i] < 100.):
                    entry += '  %8.5f' % mu[j][i]
                elif (mu[j][i] >= 100. and mu[j][i] < 1000.):
                    entry += '  %8.4f' % mu[j][i]
                elif (mu[j][i] >= 1000. and mu[j][i] < 10000.):
                    entry += '  %8.3f' % mu[j][i]
                elif (mu[j][i] >= 10000. and mu[j][i] < 100000.):
                    entry += '  %8.2f' % mu[j][i]
                else:
                    entry += '  %8.6f' % mu[j][i]

            py_printf('RESULT', entry)

        py_printf('SEPARATOR', '')

    elif isinstance(tally, Tally):

        num_bins = tallies.getNumBins()
        edges = getTallyEdges(tallies)
        mu = getTallyBatchMu(tallies)

        py_printf('HEADER', 'Batch Statistics for ' + header + ' ' + types)
        py_printf('SEPARATOR', '')

        title = ' '.center(7) + 'Energy Band' + ' '.center(9) + '    Mu    '
        py_printf('RESULT', title)

        for i in range(num_bins):
            entry = '[ %7.2f - %7.2f eV ]:' % (edges[i], edges[i+1]) 

            if (mu[i] < 1E-2):
                entry += '  ' + '%8.2E' % mu[i]
            elif (mu[i] > 10. and mu[i] < 100.):
                entry += '  %8.5f' % mu[i]
            elif (mu[i] > 100. and mu[i] < 1000.):
                entry += '  %8.4f' % mu[i]
            elif (mu[i] > 1000. and mu[i] < 10000.):
                entry += '  %8.3f' % mu[i]
            elif (mu[i] > 10000. and mu[i] < 100000.):
                entry += '  %8.2f' % mu[i]
            else:
                entry += '  ' + '%8.6f' % mu[i]

            py_printf('RESULT', entry)

        py_printf('SEPARATOR', '')

    else:
        py_printf('ERROR', 'Unable to print %s since the input is of type %s', \
                                                    types, str(type(tallies)))



###############################################################################
####################  Tally Data Processing Routines  #########################
###############################################################################

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



###############################################################################
###########################  Resonance Integrals  #############################
###############################################################################

class RIEff(object):

    def __init__(self, tally1, tally2, name=''):

        self._name=''
        self._RI=None
        self._flux=None
        self._rate=None
        self._num_RIs=0

        self.computeRIs(tally1, tally2)
        self.setName(name)


    def computeRIs(self, tally1, tally2):

        # Check that input parameters are correct
        if not isinstance(tally1, Tally):
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                   + 'integral given input of type %s', str(type(tally1)))
        if not isinstance(tally2, Tally):
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                    + 'integral given input of type %s', str(type(tally2)))
        if not tally1.hasComputedBatchStatistics():
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                    + 'integral from tally %s since it has not yet computed' \
                    + 'batch statistics', tally1.getTallyName())
        if not tally2.hasComputedBatchStatistics():
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                    + 'integral from tally %s since it has not yet computed' \
                    + 'batch statistics', tally2.getTallyName())
        if (tally1.getTallyType() is not FLUX) and \
                                        (tally2.getTallyType() is not FLUX):
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                   + 'integral since neither tally input is of FLUX tally type')

        rate_types = [CAPTURE_RATE, ELASTIC_RATE, FISSION_RATE, ABSORPTION_RATE]

        if not (tally1.getTallyType() in rate_types) and not \
                                        (tally2.getTallyType() in rate_types):
            py_printf('ERROR', 'Unable to create an effective resonance ' + \
                   'integral since neither tally input is of RATE tally type')

        if (tally1.getNumBins() is not tally2.getNumBins()):
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                       + 'integral since tally %s has %d bins while tally %s ' \
                       + 'has %d bins', tally1.getTallyName(), \
                        tally1.getNumBins(), tally2.getTallyName(), \
                        tally2.getNumBins())

        if not (getTallyEdges(tally1)==getTallyEdges(tally2)).all():
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                       + 'integral since tally %s has different bin edges ' \
                       + 'than tally %s', tally1.getTallyName(), \
                        tally2.getTallyName())


        # If we have passed all checks, then create effective RI
        self._flux = None
        self._rate = None

        if (tally1.getTallyType() is FLUX):
            self._flux = tally1
            self._rate = tally2
        else:
            self._flux = tally2
            self._rate = tally1

        self._num_RIs = self._flux.getNumBins()

        # Compute the resonance integral using tally division and 
        # multiplication operations - this allows us to compute all
        # statistics and uncertainties for the RIs, and to use all
        # of the class methods provided for a Tally class since 
        # self._RI is now a DERIVED type tally object
        edges = getTallyEdges(self._flux)
        log_array = np.log(edges[1:] / edges[0:edges.size-1])
        self._RI = (self._rate / self._flux).multiplyDoubles(log_array)


    def setName(self, name=''):
        if name is '':
            self._name = 'Effective RI'
        else:
            self._name = name

        self._RI.setTallyName(self._name)


    def getName(self):
        return self.name

    
    def getNumIntegrals(self):
        return self._RI.getNumBins()


    def getIntegrals(self):
        return getTallyBatchMu(self._RI)


    def getVariances(self):
        return getTallyBatchVariances(self._RI)


    def getStandardDeviation(self):
        return getTallyBatchStdDev(self._RI)


    def getRelativeError(self):
        return getTallyRelErr(self._RI)


    def getEnergyBandsCenters(self):
        return getTallyCenters(self._RI)


    def getEnergyBands(self):
        return getTallyEdges(self._RI)


    def getRITally(self):
        return self._RI

                
    def printRI(self, uncertainties=False):
        self._RI.printTallies(uncertainties)                


    def outputRItoFile(self, filename=''):
        self._RI.outputBatchStatistics(filename)



class RITrue(object):

    def __init__(self, isotope, bands, reaction='capture', name=''):

        if not isinstance(isotope, Isotope):
            py_printf('ERROR', 'Unable to create a true resonance integral' \
                        ' given an isotope of type %s', str(type(isotope)))

        self._name = ''
        self._isotope = isotope
        self._reaction = None
        self._num_RIs = 0
        self._num_RIs = 0
        self._RIs = None

        self.computeRIs(bands, reaction)
        self.setName(name)


    def setName(self, name=''):
        if name is '':
            self._name = 'True RI'
        else:
            self._name = name


    def computeRIs(self, bands, reaction):

        # Check that input parameters are correct
        if not isinstance(bands, list) and not isinstance(bands, np.ndarray):
            py_printf('ERROR', 'Unable to create a true resonance ' \
               + 'integral given energy bands of type %s', str(type(tally1)))

        if not isinstance(reaction, str):
            py_printf('ERROR', 'Unable to create a true resonance integral ' \
                    + 'given a reaction of type %s', str(type(reaction)))

        reaction = reaction.lower()
        rate_types = ['capture', 'elastic', 'fission', 'absorption']

        if not (reaction in rate_types):
            py_printf('ERROR', 'Unable to compute a true resonance integral ' \
                + 'for unsupported reaction rate type %s', str(type(reaction)))

        # If we have passed all checks, then create true RI
        self._energy_bands = bands
        self._reaction = reaction
        self._num_RIs = len(self._energy_bands)-1

        # Retrieve xs and energies directly from isotope
        num_energies = self._isotope.getNumXSEnergies(reaction)
        energies = self._isotope.retrieveXSEnergies(num_energies, reaction)
        xs = self._isotope.retrieveXS(num_energies, reaction)
                
        # Interpolate into the xs
        self._RIs = np.zeros(self._energy_bands.size-1)

        for i in range(self._num_RIs):
            interp_energies = np.logspace(np.log10(self._energy_bands[i]), \
                                    np.log10(self._energy_bands[i+1]), 10000)
            interp_xs = np.interp(interp_energies, energies, xs)
            self._RIs[i] = integrate.trapz(interp_xs * (1./interp_energies), \
                                                            interp_energies)


    def getName(self):
        return self.name

    
    def getNumIntegrals(self):
        return self._num_RIs


    def getIntegrals(self):
        return self._RIs


    def getEnergyBandsCenters(self):
        return getTallyCenters(self._RI)


    def getEnergyBands(self):
        return self._energy_bands

                
    def printRIs(self):

        py_printf('HEADER', 'True Resonance Integral: %s', self._name)
        py_printf('SEPARATOR', '')

        title = ' '.center(7) + 'Energy Band' + ' '.center(9) + '    RI    '
        py_printf('RESULT', title)

        for i in range(self._num_RIs):
            entry = '[ %7.2f - %7.2f eV ]:' % (self._energy_bands[i], \
                                                    self._energy_bands[i+1]) 

            if (self._RIs[i] < 1E-2):
                entry += '  ' + '%8.2E' % self._RIs[i]
            else:
                entry += '  ' + '%8.6f' % self._RIs[i]

            py_printf('RESULT', entry)

        py_printf('SEPARATOR', '')



def printRIs(RIs, header=''):
    if isinstance(RIs, list):
        if isinstance(RIs[0], RIEff):
            printTallies(RIs, header, types='Effective RIs')
        elif isinstance(RIs[0], RITrue):
            printTallies(RIs, header, types='True RIs')
        else:
            py_printf('ERROR', 'Unable to print RIs since input contains' + \
                    'contains elements of type %s', str(type(RIs)))
    elif isinstance(RIs, RIEff):
        printTallies(RIs, header, types='Effective RIs')
    elif isinstance(RIs, RITrue):
        printTallies(RIs, header, types='True RIs')
    else:
        py_printf('ERROR', 'Unable to print RIs since input is of type %s', \
                                                     str(type(RIs)))

# TODO: Plotting functions


###############################################################################
############################  Group Cross-Sections  ###########################
###############################################################################

class GroupXS(object):
    
    def __init__(self, tally1, tally2, name=''):
        
        self._name=''
        self._xs=None
        self._flux=None
        self._rate=None
        self._num_xs=0

        self.computeGroupXS(tally1, tally2)
        self.setName(name)


    def computeGroupXS(self, tally1, tally2):

        # Check that input parameters are correct
        if not isinstance(tally1, Tally):
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                   + 'given input of type %s', str(type(tally1)))
        if not isinstance(tally2, Tally):
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                    + 'given input of type %s', str(type(tally2)))
        if not tally1.hasComputedBatchStatistics():
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                    + 'from tally %s since it has not yet computed' \
                    + 'batch statistics', tally1.getTallyName())
        if not tally2.hasComputedBatchStatistics():
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                    + 'from tally %s since it has not yet computed' \
                    + 'batch statistics', tally2.getTallyName())
        if (tally1.getTallyType() is not FLUX) and \
                                        (tally2.getTallyType() is not FLUX):
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                       + 'since neither tally input is of FLUX type')

        rate_types = [CAPTURE_RATE, ELASTIC_RATE, \
                        FISSION_RATE, ABSORPTION_RATE, \
                        DIFFUSION_RATE, TRANSPORT_RATE, COLLISION_RATE ]

        if not (tally1.getTallyType() in rate_types) and not \
                                        (tally2.getTallyType() in rate_types):
            py_printf('ERROR', 'Unable to create a group cross-section ' + \
                       'since neither tally input is of RATE tally type')

        if (tally1.getNumBins() is not tally2.getNumBins()):
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                       + 'since tally %s has %d bins while tally %s ' \
                       + 'has %d bins', tally1.getTallyName(), \
                        tally1.getNumBins(), tally2.getTallyName(), \
                        tally2.getNumBins())

        if not (getTallyEdges(tally1)==getTallyEdges(tally2)).all():
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                       + 'since tally %s has different bin edges ' \
                       + 'than tally %s', tally1.getTallyName(), \
                        tally2.getTallyName())


        # If we have passed all checks, then create group XS
        self._flux = None
        self._rate = None

        if (tally1.getTallyType() is FLUX):
            self._flux = tally1
            self._rate = tally2
        else:
            self._flux = tally2
            self._rate = tally1

        self._num_xs = self._flux.getNumBins()

        # Compute the group cross-section using tally division and 
        # multiplication operations - this allows us to compute all
        # statistics and uncertainties for the xs's, and to use all
        # of the class methods provided for a Tally class since 
        # self._xs is now a DERIVED type tally object
        self._xs = (self._rate / self._flux)


    def setName(self, name=''):
        if name is '':
            if self._rate.getTallyType() is ELASTIC_RATE:
                self._name = 'Elastic Group XS'
            elif self._rate.getTallyType() is CAPTURE_RATE:
                self._name = 'Capture Group XS'
            elif self._rate.getTallyType() is FISSION_RATE:
                self._name = 'Fission Group XS'
            elif self._rate.getTallyType() is ABSORPTION_RATE:
                self._name = 'Absorption Group XS'
            elif self._rate.getTallyType() is TRANSPORT_RATE:
                self._name = 'Transport Group XS'
            elif self._rate.getTallyType() is DIFFUSION_RATE:
                self._name = 'Diffusion Coeff.'
            elif self._rate.getTallyType() is COLLISION_RATE:
                self._name = 'Total Group XS'
        else:
            self._name = name

        self._xs.setTallyName(self._name)


    def getName(self):
        return self.name

    
    def getNumXS(self):
        return self._xs.getNumBins()


    def getXS(self):
        return getTallyBatchMu(self._xs)


    def getVariances(self):
        return getTallyBatchVariances(self._xs)


    def getStandardDeviation(self):
        return getTallyBatchStdDev(self._xs)


    def getRelativeError(self):
        return getTallyRelErr(self._xs)


    def getEnergyBandsCenters(self):
        return getTallyCenters(self._xs)


    def getEnergyBands(self):
        return getTallyEdges(self._xs)


    def getRITally(self):
        return self._xs

                
    def printXS(self, uncertainties=False):
        self._xs.printTallies(uncertainties)


    def outputXStoFile(self, filename=''):
        self._xs.outputBatchStatistics(filename)
