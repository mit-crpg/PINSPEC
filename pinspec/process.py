##
# @file process.py
# @package pinspec.process
# @brief The process module provides utility functions to retrieve data
#        from PINSPEC's C++ classes, in particular, tally data. In addition,
#        the process module includes routines to use tallies to compute
#        resonance integrals and group cross-sections and to print the results
#        to the screen.
# 
# @author William Boyd (wboyd@mit.edu)
# @author Samuel Shaner (shaner@mit.edu)
# @date March 25, 2013

import numpy as np
import matplotlib.pyplot as plt
import pinspec
from log import *
import scipy.integrate as integrate


##
# @brief Returns an array of the center values for a tally's bins.
# @details A wrapper function to make it easier to access a tally's bin center
#          data array through SWIG. This would be invoked a PINSPEC Python
#          input file as follows:
#
# @code
#          bin_center_array = pinspec.process.getTallyCenters(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally bin centers
def getTallyCenters(tally):

    if not isinstance(tally, pinspec.Tally):
        py_printf('WARNING', 'Unable to get tally centers from input of type' \
                      + str(type(tally)))
    else:
        num_bins = tally.getNumBins()
        centers = tally.retrieveTallyCenters(num_bins)
        return centers


##
# @brief Returns an array of the bin edge values for a tally's bins.
# @details A wrapper function to make it easier to access a tally's bin edges
#          data array through SWIG. This would be invoked a PINSPEC Python
#          input file as follows:
#
# @code
#          bin_edges_array = pinspec.process.getTallyEdges(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally bin edges
def getTallyEdges(tally):

    if not isinstance(tally, pinspec.Tally):
        py_printf('WARNING', 'Unable to get tally edges from input of type' \
                                                            + str(type(tally)))
    else:
        num_bins = tally.getNumBins()
        edges = tally.retrieveTallyEdges(num_bins+1)
        return edges


##
# @brief Returns an array of the batch averages for the tally's bins.
# @details A wrapper function to make it easier to access a tally's tallies
#          data array through SWIG. This would be invoked a PINSPEC Python
#          input file as follows:
#
# @code
#          tally_averages = pinspec.process.getTallyBatchMu(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally batch averages
def getTallyBatchMu(tally):

    if not isinstance(tally, pinspec.Tally):
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
    

##
# @brief Returns an array of the batch variances for a tally's bins.
# @details A wrapper function to make it easier to access a tally's variances
#          data array through SWIG. This would be invoked a PINSPEC Python
#          input file as follows:
#
# @code
#          tally_variances = pinspec.process.getTallyBatchVariances(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally batch variances
def getTallyBatchVariances(tally):

    if not isinstance(tally, pinspec.Tally):
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


##
# @brief Returns an array of the batch standard deviations for a tally's bins.
# @details A wrapper function to make it easier to access a tally's standard
#          deviations data array through SWIG. This would be invoked a PINSPEC 
#          Python input file as follows:
#
# @code
#          tally_std_dev = pinspec.process.getTallyBatchStdDev(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally batch standard deviations
def getTallyBatchStdDev(tally):

    if not isinstance(tally, pinspec.Tally):
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


##
# @brief Returns an array of the batch relative errors for a tally's bins.
# @details A wrapper function to make it easier to access a tally's relative
#          errors data array through SWIG. This would be invoked a PINSPEC 
#          Python input file as follows:
#
# @code
#          tally_rel_err = pinspec.process.getTallyBatchRelErr(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally batch relative errors
def getTallyBatchRelErr(tally):

    if not isinstance(tally, pinspec.Tally):
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


##
# @brief Returns a 2D array of the tally's batch statistics.
# @details This is a wrapper function to make it easier to access a tally's
#          batch-based statistical data through SWIG. The array returned 
#          contains the tally bin centers, averages, variances, standard 
#          deviations, and relative errors, in that order. This would be 
#          invoked a PINSPEC Python input file as follows:
#
# @code
#          tally_data = pinspec.process.getTallyBatchStatistics(tally)
#          bin_centers = tally_data[0][:]
#          tally_avareages = tally_data[1][:]
#          tally_variances = tally_data[2][:]
#          tally_std_dev = tally_data[3][:]
#          tally_rel_err = tally_data[4][:]
# @endcode
#
# @param tally the tally of interest
# @return a 2D numpy array with the tally batch statistical data
def getTallyBatchStatistics(tally):

    if not isinstance(tally, pinspec.Tally):
        py_printf('WARNING', 'Unable to get tally statistics from input of ' \
                                                + ' type' + str(type(tally)))
    if not tally.hasComputedBatchStatistics():
        py_printf('WARNING', 'Unable to get tally statistics for tally %s'
                        ' since it has not yet computed batch statistics', \
                       + tally.getTallyName())
    else:
        centers = getTallyCenters(tally)
        mu = getTallyBatchMu(tally)
        variances = getTallyBatchVariances(tally)
        std_dev = getTallyBatchStdDev(tally)
        rel_err = getTallyBatchRelErr(tally)

        # Create a 2D array of all arrays        
        statistics = np.array([centers, mu, variances, std_dev, rel_err])
        return statistics


##
# @brief Prints formatted tally batch-averaged data to the screen as a table.
# @details Prints a formatted table of tally data to the screen and can be used
#          for a single tally or for a list of tallies. Since RIEff objects are 
#          simply Python wrappers for an underlying DERIVED type tally, this 
#          method can also be used to print lists of RIEff objects. It also 
#          works for RITrue objects, though RITrue objects are not stored as
#          tallies. A user may invoke this function from a PINSPEC Python input 
#          file as follows:
#
# @code
#          printTallies([flux1, flux2, flux3], header='Flux Tallies')
# @endcode
#
# @param tallies a list of the tallies to print to the screen
# @param header an optional string to prepend to the title of table.
# @param types an optional string of the tally types for the table title
#
def printTallies(tallies, header='', types='Tallies'):

    if type(tallies) is list:

        # Check that all elements are Tally ojbects
        for tally in tallies:

            if isinstance(tally, RIEff):
                for i in range(len(tallies)):
                    tallies[i] = tallies[i].getRITally()

            elif not isinstance(tally, pinspec.Tally) and not isinstance(tally, RITrue):
                py_printf('ERROR', 'Unable to print %s since input contains' + \
                        'contains elements of type %s', types, str(type(tally)))

        if isinstance(tallies[0], pinspec.Tally):

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

    elif isinstance(tally, pinspec.Tally):

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


'''
Tally Data Processing Routines
'''

##
# @brief Computes the mean number of collisions for a neutron before absorption.
# @param coll_rate a COLLISION_RATE tally
# @param num_neutrons the number of neutrons per batch
# @return the mean number of collisions per neutron before absorption
#
def computeMeanNumCollisions(coll_rate, num_neutrons):

    coll_rate.computeScaledBatchStatistics(num_neutrons)

    num_bins = coll_rate.getNumBins()
    coll_rate_mu = coll_rate.retrieveTallyMu(num_bins)
    mean_rate = 0.0

    for i in range(num_bins):
        mean_rate += coll_rate_mu[i]

    return mean_rate


##
# @brief Computes the mean neutron lifetime (seconds) before absorption.
# @param coll_times INTERCOLLISION_TIME tally
# @param num_neutrons the number of neutrons per batch
# @return the mean neutron lifetime (seconds) before absoprtion
def computeMeanNeutronLifetime(coll_times, num_neutrons):

    coll_times.computeScaledBatchStatistics(num_neutrons)

    num_bins = coll_times.getNumBins()
    coll_times_mu = coll_times.retrieveTallyMu(num_bins)
    mean_time = 0.0

    for i in range(num_bins):
        mean_time += coll_times_mu[i]

    return mean_time


##
# @class RIEff process.py "pinspec/process.py"
# @brief An effective resonance integral.
# @details This class represents an effective resonance integral computed
#          by a PINSPEC monte carlo simulation. Two tallies - one reaction
#          rate and one flux - are required to form a RIEff class object.
# 
class RIEff(object):

    ##
    # @brief RIEff constructor.
    # @param self the RIEff object pointer
    # @param tally1 one of the tallies needed to compute an \f$ RI_{eff} \f$
    # @param tally2 one of the tallies needed to compute an \f$ RI_{eff} \f$
    # @param name an optional string for the name of the resonance integral
    def __init__(self, tally1, tally2, name=''):

        ## The name of the effective resonance integral
        self._name=''
        ## The DERIVED tally type for the effective resonance integral data
        self._RI=None
        ## The pinspec.FLUX tally used to compute the effective resonance integral
        self._flux=None
        ## The REACTION_RATE tally used to compute the effective 
        #  resonance integral 
        self._rate=None
        ## The number of resonance integral energy bands
        self._num_RIs=0

        self.computeRIs(tally1, tally2)
        self.setName(name)


    ##
    # @brief Computes the resonance integrals from the two tallies.
    # @details This method checks that one of the tallies input is a reaction
    #          rate and the other is a flux, with the same bin edges. If the
    #          batch statistics for both tallies have been computed, it computes
    #          the resonance integral as follows: 
    #
    # @code
    #          edges = getTallyEdges(tally1)
    #          log_array = numpy.log(edges[1:] / edges[0:edges.size-1])
    #          RI_eff = (reaction_rate / flux).multiplyDoubles(log_array)
    # @endcode
    #
    # @param self the RIEff object pointer
    # @param tally1 one of the tallies needed to compute an \f$ RI_{eff} \f$
    # @param tally2 one of the tallies needed to compute an \f$ RI_{eff} \f$
    def computeRIs(self, tally1, tally2):

        # Check that input parameters are correct
        if not isinstance(tally1, pinspec.Tally):
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                   + 'integral given input of type %s', str(type(tally1)))
        if not isinstance(tally2, pinspec.Tally):
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
        if (tally1.getTallyType() is not pinspec.FLUX) and \
                                        (tally2.getTallyType() is not pinspec.FLUX):
            py_printf('ERROR', 'Unable to create an effective resonance ' \
                   + 'integral since neither tally input is of pinspec.FLUX tally type')

        rate_types = [pinspec.CAPTURE_RATE, \
                      pinspec.ELASTIC_RATE, \
                      pinspec.FISSION_RATE, \
                      pinspec.ABSORPTION_RATE]

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

        if (tally1.getTallyType() is pinspec.FLUX):
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


    ##
    # @brief Sets the name of this effective resonance integral.
    # @details This is useful when one wishes to print the resonance integral
    #          values to the screen or a file since it will be identifiable
    #          by the user-defined name.
    # @param self the RIEff object pointer
    # @param name the name of the RIEff object
    def setName(self, name=''):
        if name is '':
            self._name = 'Effective RI'
        else:
            self._name = name

        self._RI.setTallyName(self._name)


    ##
    # @brief Returns the name of the RIEff object.
    # @details Returns an empty string if no name has been specified by the user.
    # @param self the RIEff object pointer
    # @return a string with the name of the RIEff
    def getName(self):
        return self.name

    ##
    # @brief Returns the number of resonance integral energy bands.
    # @details The number of integrals is equivalent to the number of 
    #          tally bins for the reaction rate and flux tallies.
    # @param self the RIEff object pointer
    # @return the number of resonance integrals
    def getNumIntegrals(self):
        return self._RI.getNumBins()


    ##
    # @brief Returns an array of the resonance integral batch averaged values
    #        for each tally bin.
    # @param self the RIEff object pointer
    # @return a numpy array of the resonance integrals
    def getIntegrals(self):
        return getTallyBatchMu(self._RI)


    ##
    # @brief Returns an array of the resonance integral batch variances for each
    #        tally bin.
    # @param self the RIEff object pointer
    # @return a numpy array of the resonance integral variances
    def getVariances(self):
        return getTallyBatchVariances(self._RI)


    ##
    # @brief Returns an array of the resonance integral batch standard deviations
    #        for each tally bin.
    # @param self the RIEff object pointer
    # @return a numpy array of the resonance integral standard deviations
    def getStandardDeviation(self):
        return getTallyBatchStdDev(self._RI)

 
    ##
    # @brief Returns an array of the resonance integral batch relative errors
    #        for each tally bin.
    # @param self the RIEff object pointer
    # @return a numpy array of the resonance integral relative errors
    def getRelativeError(self):
        return getTallyBatchRelErr(self._RI)


    ##
    # @brief Returns an array of the resonance integral energy band centers
    #        for each tally bin.
    # @param self the RIEff object pointer
    # @return a numpy array of the resonance integral energy band centers
    def getEnergyBandsCenters(self):
        return getTallyCenters(self._RI)


    ##
    # @brief Returns an array of the resonance integral energy bands.
    # @details The energy bands are the tally bin edges for the reaction rate
    #          and flux tally forming this effective resonance integral.
    # @param self the RIEff object pointer
    # @return a numpy array of the resonance integral energy bands (eV)
    def getEnergyBands(self):
        return getTallyEdges(self._RI)


    ##
    # @brief Returns a reference to this RIEff object.
    # @param self RIEff the object pointer
    # @return a reference to the RIEff object
    def getRITally(self):
        return self._RI

                
    ##
    # @brief Prints a formatted table of the effective resonance integrals to 
    #        the screen.
    # @details The resonance integrals and their uncertainties (optional) will
    #          be printed as a formatted table to the screen.
    # @param self the RIEff object pointer
    # @param uncertainties whether or not to print tally statistics (default is
    #        false)
    # @return a reference to the RIEff object
    def printRI(self, uncertainties=False):
        self._RI.printTallies(uncertainties)                


    ##
    # @brief Prints the resonance integral data and batch statistics to a file.
    # @details Since the effective resonance integral is stored as a DERIVED 
    #          tally type, this method prints the resonance integral data to a
    #          file using the Tally::outputBatchStatistics() method. An
    #          auto-generated filename will be created with the format
    #          'tally-#.data' where # is an auto-incremented integer for each
    #          tally output data file created.
    # @param self the RIEff object pointer
    # @param filename An optional filename for the output file
    def outputRItoFile(self, filename=''):
        self._RI.outputBatchStatistics(filename)



##
# @class RITrue process.py "pinspec/process.py"
# @brief A true resonance integral.
# @details This class represents the true resonance integral computed
#          from a \f$ \frac{1}{E} \f$ spectrum.
class RITrue(object):

    ##
    # @brief RITrue constructor.
    # @param self the RITrue object pointer
    # @param isotope a pointer to the isotope of interest
    # @param bands an array of the energy bands for each resonance integral
    # @param reaction an optional argument string with the reaction rate type
    # @param name an optional argument string for the resonance integral name
    def __init__(self, isotope, bands, reaction='capture', name=''):

        if not isinstance(isotope, pinspec.Isotope):
            py_printf('ERROR', 'Unable to create a true resonance integral' \
                        ' given an isotope of type %s', str(type(isotope)))

        ## The name of the effective resonance integral
        self._name = ''
        ## The isotope for this resonance integral
        self._isotope = isotope
        ## The reaction rate for this resonance integral (ie, 'capture')
        self._reaction = None
        ## The number of resonance integrals
        self._num_RIs = 0
        ## The numpy array of resonance integral values
        self._RIs = None

        self.computeRIs(bands, reaction)
        self.setName(name)


    ##
    # @brief Sets the name of the true resonance integral.
    # @details This is useful when one wishes to print the resonance integral
    #          values to the screen or a file since it will be identifiable
    #          by the user-defined name.
    # @param self the RITrue object pointer
    # @param name the name of the RITrue object
    def setName(self, name=''):
        if name is '':
            self._name = 'True RI'
        else:
            self._name = name


    ##
    # @brief Computes the true resonance integrals for an isotopic reaction 
    #        rate.
    # @details This method computes an infinite dilute resonance integral for
    #          each of the user-defined energy bands for some isotopic reaction 
    #          rate using numerical integration as follows:
    #
    #      \f$ \int\limits_{E_{i-1}}^{E_i} \sigma_{j}(E) \frac{1}{E} 
    #      \mathrm{d}E \f$
    #          
    # @param self the RITrue object pointer
    # @param bands an array of the energy bands for each resonance integral
    # @param reaction an optional argument string with the reaction rate type
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
        ## The numpy array of energy band values
        self._energy_bands = bands
        ## The reaction rate type (ie, 'capture')
        self._reaction = reaction
        ## The number of resonance integrals
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


    ##
    # @brief Returns the name of the RITrue object.
    # @details Returns an empty string if no name has been specified by the 
    #          user.
    # @param self the RITrue object pointer
    # @return a string with the name of the RITrue
    def getName(self):
        return self.name

    
    ##
    # @brief Returns the number of resonance integral energy bands.
    # @param self the RITrue object pointer
    # @return the number of resonance integrals
    def getNumIntegrals(self):
        return self._num_RIs


    ##
    # @brief Returns an array of the resonance integrals for each energy band.
    # @param self the RITrue object pointer
    # @return a numpy array of the resonance integrals
    def getIntegrals(self):
        return self._RIs


    ##
    # @brief Returns an array of the centers of each energy band.
    # @param self the RITrue object pointer
    # @return a numpy array of the energy band centers (eV)
    def getEnergyBandsCenters(self):
        return getTallyCenters(self._RI)


    ##
    # @brief Returns an array of the energy bands.
    # @param self the RITrue object pointer
    # @return a numpy array of the values defining the energy bands (eV)
    def getEnergyBands(self):
        return self._energy_bands

                
    ##
    # @brief Prints a formatted table of the true resonance integrals to the
    #        screen.
    # @param self the RITrue object pointer
    # @return a reference to the RITrue object
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



##
# @brief Prints a formatted table for an array of true and/or effectrive 
#        resonance integrals to the screen.
# @param RIs a list of resonance integrals (RIEff or RITrue objects)
# @param header an optional argument string for the table title
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


##
# @class GroupXS process.py "pinspec/process.py"
# @brief A multi-group cross-section for a certain reaction rate.
# @details This class represents a multi-group cross-section with some energy
#          bands for a certain reaction rate. By default, the group 
#          cross-section is a macroscopic cross-section \f$ (cm^{-1}) \f$.
class GroupXS(object):
    
    ##
    # @brief The group cross-section constructor.
    # @param self the GroupXS object pointer
    # @param tally1 one of the tallies needed to compute a multi-group 
    #         cross-section  
    # @param tally2 one of the tallies needed to compute a multi-group 
    #        cross-section
    # @param name an optional argument string for the name of the multi-group
    #        cross-section
    def __init__(self, tally1, tally2, name=''):
        
        ## The name of the multi-group cross-section
        self._name=''
        ## The DERIVED tally used to store the multi-group cross-sections
        self._xs=None
        ## The pinspec.FLUX tally used to compute the multi-group cross-sections
        self._flux=None
        ## The REACTION_RATE tally used to compute the multi-group 
        #  cross-sections
        self._rate=None
        ## The number of multi-group cross-sections
        self._num_xs=0

        self.computeGroupXS(tally1, tally2)
        self.setName(name)

 
    ##
    # @brief Computes the multigroup cross-sections from the two tallies.
    # @details This method checks that one of the tallies input is a reaction
    #          rate and the other is a flux, with the same bin edges. If the
    #          batch statistics for both tallies have been computed, it computes
    #          the multigroup cross-section as follows: 
    #
    # @code
    #          sigma_groups = reaction_rate / flux
    # @endcode
    #
    # @param self the GroupXS object pointer
    # @param tally1 one of the tallies needed to compute a multi-group 
    #        cross-section
    # @param tally2 one of the tallies needed to compute a multi-group
    #        cross-section
    def computeGroupXS(self, tally1, tally2):

        # Check that input parameters are correct
        if not isinstance(tally1, pinspec.Tally):
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                   + 'given input of type %s', str(type(tally1)))
        if not isinstance(tally2, pinspec.Tally):
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
        if (tally1.getTallyType() is not pinspec.FLUX) and \
           (tally2.getTallyType() is not pinspec.FLUX):
            py_printf('ERROR', 'Unable to create a group cross-section ' \
                       + 'since neither tally input is of pinspec.FLUX type')

        rate_types = [pinspec.CAPTURE_RATE, \
                      pinspec.ELASTIC_RATE, \
                      pinspec.FISSION_RATE, \
                      pinspec.ABSORPTION_RATE, \
                      pinspec.DIFFUSION_RATE, \
                      pinspec.TRANSPORT_RATE, \
                      pinspec.COLLISION_RATE]

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

        if (tally1.getTallyType() is pinspec.FLUX):
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


    ##
    # @brief Sets the name of this group cross-section.
    # @details This is useful when one wishes to print the multi-group 
    #          cross-section values to the screen or a file since it will be 
    #          identifiable by the user-defined name.
    # @param self the GroupXS object pointer
    # @param name the name of the GroupXS object
    def setName(self, name=''):
        if name is '':
            if self._rate.getTallyType() is pinspec.ELASTIC_RATE:
                self._name = 'Elastic Group XS'
            elif self._rate.getTallyType() is pinspec.CAPTURE_RATE:
                self._name = 'Capture Group XS'
            elif self._rate.getTallyType() is pinspec.FISSION_RATE:
                self._name = 'Fission Group XS'
            elif self._rate.getTallyType() is pinspec.ABSORPTION_RATE:
                self._name = 'Absorption Group XS'
            elif self._rate.getTallyType() is pinspec.TRANSPORT_RATE:
                self._name = 'Transport Group XS'
            elif self._rate.getTallyType() is pinspec.DIFFUSION_RATE:
                self._name = 'Diffusion Coeff.'
            elif self._rate.getTallyType() is pinspec.COLLISION_RATE:
                self._name = 'Total Group XS'
        else:
            self._name = name

        self._xs.setTallyName(self._name)


    ##
    # @brief Returns the name of the GroupXS object.
    # @details Returns an empty string if no name has been specified by the 
    #          user.
    # @param self the GroupXS object pointer
    # @return a string with the name of the GroupXS
    def getName(self):
        return self.name


    ##
    # @brief Returns the number of multi-group cross-section values
    # @param self the GroupXS object pointer
    # @return the number of multi-group cross-section values
    def getNumXS(self):
        return self._xs.getNumBins()


    ##
    # @brief Retrurns an array of the batch-averaged multi-group cross-sections.
    # @param self the GroupXS object pointer
    # @return a numpy array of the multi-group cross-sections
    def getXS(self):
        return getTallyBatchMu(self._xs)


    ##
    # @brief Returns an array of the multi-group cross-section variances.
    # @param self the GroupXS object pointer
    # @return a numpy array of the multi-group cross-section variances
    def getVariances(self):
        return getTallyBatchVariances(self._xs)


    ##
    # @brief Returns an array of the multi-group cross-section standard 
    #        deviations.
    # @param self the GroupXS object pointer
    # @return a numpy array of the multi-group cross-section standard deviations
    def getStandardDeviation(self):
        return getTallyBatchStdDev(self._xs)


    ##
    # @brief Returns an array of the multi-group cross-section relative errors.
    # @param self the GroupXS object pointer
    # @return a numpy array of the multi-group cross-section relative errors
    def getRelativeError(self):
        return getTallyBatchRelErr(self._xs)


    ## 
    # @brief Returns an array of the multi-group cross-section energy 
    #        band centers.
    # @param self the GroupXS object pointer
    # @return a numpy array of the multi-group cross-section energy band centers
    def getEnergyBandsCenters(self):
        return getTallyCenters(self._xs)


    ## 
    # @brief Returns an array of the multi-group cross-section energy
    #        band values.
    # @param self the GroupXS object pointer
    # @return a numpy array of the multi-group cross-section energy band values
    def getEnergyBands(self):
        return getTallyEdges(self._xs)


    ##
    # @brief Returns a reference to this GroupXS object.
    # @param self GroupXS the object pointer
    # @return a reference to the GroupXS object
    def getGroupXS(self):
        return self._xs

                
    ##
    # @brief Prints a formatted table of the multi-group cross-sections to
    #        the screen.
    # @details The multi-group cross-sections and their uncertainties 
    #          (optional) will be printed as a formatted table to the screen.
    # @param self the GroupXS object pointer
    # @param uncertainties whether or not to print tally statistics 
    #        (default is false)
    # @return a reference to the GroupXS object
    def printXS(self, uncertainties=False):
        self._xs.printTallies(uncertainties)


    ##
    # @brief Prints the multi-group cross-section data and batch statistics 
    #        to a file.
    # @details Since the multi-group cross-sections are stored as a DERIVED 
    #          tally type, this method prints the cross-section data to a
    #          file using the Tally::outputBatchStatistics() method. An
    #          auto-generated filename will be created with the format
    #          'tally-#.data' where # is an auto-incremented integer for each
    #          tally output data file created.
    # @param self the GroupXS object pointer
    # @param filename An optional filename for the output file
    def outputXStoFile(self, filename=''):
        self._xs.outputBatchStatistics(filename)
