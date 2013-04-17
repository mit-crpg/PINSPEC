##
# @file slbw.py
# @package pinspec.slbw
# @brief The slbw module provides utility functions to generate resonant 
#        cross-sections using the Single-Level Breit-Wigner formalism. This
#        is most helpful for PINSPEC applications which are temperature-
#        dependent and require doppler broadening effects to be treated for
#        resonant absorbers.
# @author Jessica Hunter
# @date April 17, 2013

from pinspec import *
from log import *
import pinspec
import matplotlib.pyplot as plt
import numpy
import os
import shutil
import subprocess
import math
import scipy.special as spec

##
# @brief Function to create resonant capture and scatter cross-sections.
# @details Generates the cross-section at some temperature using resonance
#          data and the Single-Level Breit-Wigner formalism. Function broadens
#          a specific number of resonances from the ENDFB-VII library, then 
#          creates identical resonances with even spacing to a specified
#          energy limit, and then a specific flat cross section for the 
#          remaining energy. When generating elastic cross sections, the
#          flat cross section is set to the potential scattering value.
# @code
# temp=1200
# slbw.SLBWXS(u238, temp, 'capture')
# @endcode
# @param isotope the isotope of interest
# @param temp the temperature in degrees Kelvin
# @param xs_type an optional argument string for the cross-section type
# @param number_of_pos_res the number of positive resonance to use with
#        \f$ E_0 > 0 \f$
# @param energy_min the minimum energy at which to generate the cross-section (eV)
# @param energy_max the maximum energy of the cross section (eV)
# @param upper_energy_limit_identical_res the upper energy limit of the equally spaced, identical resonances section
# @param energy_bin_width the width of the energy bin for the cross section
# @param resonance_spacing_identical_res the spacing between resonance peaks for the identical resonances section
# @param lower_bound_identical_res the lower energy bound for the identical resonances section
# @param gamma_gamma_identical_res the \f$\Gamma_{\gamma}\f$ value for the identical resonances section
# @param flat_xs a flat cross section used for the capture cross section above the identical resonances region
def buildSLBWXS(isotope,temp,xs_type='capture',number_of_pos_res=14, \
           energy_min=1e-5, energy_max=20E6, upper_energy_limit_identical_res=1000.0, energy_bin_width = 0.075, resonance_spacing_identical_res = 25, \
           lower_bound_identical_res = 300, gamma_gamma_identical_res = 0.023, flat_xs = 0.1):

    #-------------------------------------------------
    #---Construct resonance parameter arrays from given file
    #--------------------------------------------------
    # Get file path 
    filepath = str(pinspec.getXSLibDirectory()) + isotope +'-RP.txt'
    if os.path.exists(filepath) == True:
 	py_printf('INFO', 'Loading resonance paramater file '\
 						 + isotope + '-RP.txt')
    else:
 	py_printf('WARNING', 'Unable to load resonance parameter '\
                                          + isotope + '-RP.txt')
    # Initialize desired arrays
    E0 = numpy.array([])
    GN = numpy.array([])
    GG = numpy.array([])
    with open(filepath) as restxt:
	# Parse first line for SigP
	junk, SigP, barns = restxt.readline().split(' ', 2)
	SigP = float(SigP)
	# Skip formatting lines
	next(restxt)
	next(restxt)
	for line in restxt:
 	    # See if first line is negative
	    num, restofline = line.split('.', 1)
	    num = float(num)
	    # skip negative resonances
	    if num < 0:
		next(restxt)
	    else:
		break
        # For first positive resonance, parse line back together
        num = int(num)
	num = str(num)
	firstposresonance = num + '.' + restofline
	# Parse first line properly, then append
	firstline = parse(firstposresonance)
	E0 = numpy.append(E0, firstline[0])
	GN = numpy.append(GN, firstline[1])
	GG = numpy.append(GG, firstline[2])
	# parse each of the rest of the desired resonances and append to 
        # the arrays
	for i, line in enumerate(restxt):
            if i < (number_of_pos_res - 1):
                out = parse(line)
		E0 = numpy.append(E0, out[0])
		GN = numpy.append(GN, out[1])
		GG = numpy.append(GG, out[2])
            else:
                break

    restxt.close()
    #--------------------------------------------------
    # Now we have SigP, E0, GN, and GG for first 14 (for U238) resonances, 
    # need  rest
    #--------------------------------------------------
    #resonances up to 1keV------------
    # Make up an E_0 vector with resonance spacing from 300 to 1 keV
    numbns = round((upper_energy_limit_identical_res - lower_bound_identical_res) / resonance_spacing_identical_res)
    Enot = numpy.linspace(lower_bound_identical_res, upper_energy_limit_identical_res, numbns)

    # Make up Gamma_gamma (0.023 for U-238)
    Gamma_g = numpy.empty_like(Enot)
    Gamma_g[:] = gamma_gamma_identical_res

    # Make up gamma_n for each
    # Last known resonance E_0 (291.0206 for U-238)
    E_last = E0[number_of_pos_res - 1]  
    Gamma_nlast = GN[number_of_pos_res - 1]  # Last known Gn value
    Gamma_n = numpy.empty_like(Enot)
    for k in range(len(Enot)):
	Gamma_n[k] = Gamma_nlast * ((Enot[k] / E_last) ** 0.5)
	E_last = Enot[k]

    # Append arrays together
    E0 = numpy.append(E0, Enot)
    GN = numpy.append(GN, Gamma_n)
    GG = numpy.append(GG, Gamma_g)
    #--------------------------------------------------
    # SLBW section
    #--------------------------------------------------

    # pull A value from isotope
    El, A = isotope.split('-', 1)
    A = float(A)

    # Construct Energy grid for fictitious XS
    nbins = round((upper_energy_limit_identical_res - energy_min) / energy_bin_width)
    E = numpy.logspace(numpy.log10(energy_min), numpy.log10(upper_energy_limit_identical_res), nbins)

    # Set k value
    k = 8.617e-5

    # allocate and clear arrays
    x = numpy.empty_like(E)
    sigma_g = numpy.zeros_like(E)
    sigma_n = numpy.zeros_like(E)
    sigma_n[:] = SigP
    psichi = numpy.empty_like(E)

    # arrays for each resonance (CHECK IF NUMPY MATH IS RIGHT)
    Gamma = GN + GG  # total gamma
    squig = Gamma * ((A / (4 * k * temp * E0)) ** 0.5)  # squiggle
    r = (2603911 / E0) * ((A + 1) / A)  # r
    q = (r * SigP) ** 0.5  # q

    # generate cross sections
    for i in range(len(E0)):  # for each resonance

	# create little chi or x vector for a single resonance
	for j in range(len(E)):
 	    x[j] = 2 * (E[j] - E0[i]) / Gamma[i]

	# use W code to generate psi and chi functions
	y = (((x + 1j) / 2) * squig[i])
	psichi = (math.pi * squig[i] / (2 * (math.pi ** 0.5))) * spec.wofz(y)
	psi = psichi.real
	chi = 2.0*psichi.imag

	# generate cross sections for this resonance, and sum
	sigma_g = sigma_g + (GN[i] / Gamma[i]) * (GG[i] / Gamma[i]) * \
	    ((E0[i] / E) ** 0.5) * (r[i] * psi)
	sigma_n = sigma_n + ((GN[i] / Gamma[i])**2) * (r[i] * psi + q[i] * chi)		
	
    #-----------------------------------------------------------------
    # Now that we have cross sections, we can frankenstein after 1 keV
    #----------------------------------------------------------------
    # From 1keV to 20E6 eV, use 0.1 barns for capture (automated to use inputs)
    E2 = numpy.linspace(upper_energy_limit_identical_res + energy_bin_width, energy_max, 10)
    sigma_g2 = numpy.empty_like(E2)
    sigma_g2[:] = flat_xs
    sigma_n2 = numpy.empty_like(E2)
    sigma_n2[:] = SigP

    # Frankenstein arrays together
    sigma_g = numpy.append(sigma_g, sigma_g2)
    sigma_n = numpy.append(sigma_n, sigma_n2)
    E = numpy.append(E, E2)
    EAXS = (E, sigma_g)
    EAXS = numpy.transpose(EAXS)
    ESXS = (E, sigma_n)
    ESXS = numpy.transpose(ESXS)

    #----------------------------------
    #---So if capture is specified
    #---------------------------------
    if xs_type=='capture':
   	# write output file for capture XS
	out_name = getXSLibDirectory()+El+'-'+str(int(A))+'-capture.txt'
	numpy.savetxt(out_name, EAXS, newline='\n', delimiter=',')
	# go back and add header
	f = open(out_name)
	text = f.read()
	f.close()
	# open the file again for writing
	f = open(out_name, 'w')
	f.write('Doppler Broadened SLBW fictitious capture XS at temp=' + \
			str(temp)+'K\n')
	# write the original contents
	f.write(text)
	f.close()
    #----------------------------------
    #---So if capture is specified
    #---------------------------------
    if xs_type=='scatter':
   	# write output file for scatter XS
	out_names = getXSLibDirectory()+El+'-'+str(int(A))+'-elastic.txt'
	numpy.savetxt(out_names, ESXS, newline='\n', delimiter=',')
	# go back and add header
	g = open(out_names)
	texts = g.read()
	g.close()
	# open the file again for writing
	g = open(out_names, 'w')
	g.write('Doppler Broadened SLBW fictitious resonant scattering XS ' + \
			'at temp='+str(temp)+'K\n')
	# write the original contents
	g.write(texts)
	g.close()
	#-----------------

##
# @brief Function to convert a string into a float.
# @param string the string we wish to convert
# @return the float version of the input string
def convert(string):
    n = 8
    num, exp = [string[i:i + n] for i in range(0, len(string), n)]
    num = float(num)
    exp = float(exp)
    value = num * (10 ** exp)
    return value

##
# @brief Function to parse a string and convert into resonance parameters.
# @details Function to parse a resonance parameter string into \f$E_0\f$,
#          \f$\Gamma_n\f$, and \f$\Gamma_{\gamma}\f$.
# @param string the string of interest
def parse(string):
    result = numpy.array([])
    E0, J, GN, GG, GFA, GFB = string.strip().split(' ', 5)
    # Further parse individual strings, convert to floats
    E0 = convert(E0)
    GN = convert(GN)
    GG = convert(GG)
    together = E0, GN, GG
    result = numpy.append(result, together)
    return result


##
# @brief Function to generate a potential scattering cross-section
#        for an isotope based on data in a resonance parameters file.
# @param isotope the isotope of interest
# @param energy_min the upper limit of the energy range
# @param energy_max the lower limit of the energy range
def generatePotentialScattering(isotope, energy_min=1e-5, energy_max=20E6):

    #isotope = 'U-238-ResonanceParameters.txt'
    filepath = str(pinspec.getXSLibDirectory()) + isotope + '-RP.txt'
    if os.path.exists(filepath) == True:
	py_printf('INFO', 'Loading resonance paramater file' + filepath)
    else:
	py_printf('WARNING', 'Unable to load resonance ' \
                                             'parameter file' + filepath)
    with open(filepath) as restxt:
	# Parse first line for SigP
	junk, SigP, barns = restxt.readline().split(' ', 2)
	SigP = float(SigP)
    E = numpy.logspace(numpy.log10(energy_min), numpy.log10(energy_max), 2)
    SXS=numpy.zeros_like(E)
    SXS[:]=SigP
    ESXS = (E, SXS)
    ESXS = numpy.transpose(ESXS)
    # pull A value from filename
    El, A= isotope.split('-', 1)
    A = float(A)
    # write output file for scatter XS
    out_names = getXSLibDirectory()+El+'-'+str(int(A))+'-elastic.txt'
    numpy.savetxt(out_names, ESXS, newline='\n', delimiter=',')
    # go back and add header
    g = open(out_names)
    texts = g.read()
    g.close()
    # open the file again for writing
    g = open(out_names, 'w')
    g.write('Fictitious resonant scattering XS, values are potential XS\n')
    # write the original contents
    g.write(texts)
    g.close()

##
# @brief Function to generate a plot of the SLBW generated cross-section
#        along with the ENDF-VII version of the cross-section for comparison.
# @param isotope the isotope of interest
# @param type_xs type of cross section, 'capture' or 'elastic' ('scatter' also accepted)
# @param dir the directory in which the plot will be saved
def compareXS(isotope, type_xs='capture', dir='.'):
	
    #Get fake XS from info given
    El, A = isotope.split('-', 1)
    #Find proper filename for fake XS
    if type_xs=='scatter':
		type_xs='elastic'
    path=str(getXSLibDirectory())+'/'+El+'-'+A+'-'+type_xs+'.txt'
    #Read in array for fictitious XS at 300
    EnT=numpy.array([])
    barnsT=numpy.array([])
    invEnT=numpy.array([])
    with open(path) as resT:
   	#Parse out String containing Temperature
	Junk, temp = resT.readline().split('=', 1)
	for line in resT:
 	    EnTt, barnsTt=line.split(',', 1)
	    EnTt=float(EnTt)
	    barnsTt=float(barnsTt)
	    invEnTt=1/EnTt
	    EnT=numpy.append(EnT,EnTt)
	    barnsT=numpy.append(barnsT,barnsTt)
	    invEnT=numpy.append(invEnT, invEnTt)
	
    py_printf('INFO', 'Read in Doppler Broadened XS correctly')

    #Read in array for ENDF7 XS at 300
    npath=str(getXSLibDirectory())+'/BackupXS/'+El+'-'+A+'-' + \
                                                            type_xs+'.txt'
    EndfE300=numpy.array([])
    barnsEndF300=numpy.array([])
    invEndfE300=numpy.array([])
    with open(npath) as Endf300:
	#Parse out String containing Temperature
	Junk, xssource = Endf300.readline().split(' ', 1)
	xssource=xssource.strip()
	for line in Endf300:
 	    EndfE300temp, barnsEndF300temp=line.split(',', 1)
	    EndfE300temp=float(EndfE300temp)
	    barnsEndF300temp=float(barnsEndF300temp)
	    invEndfE300temp=1/EndfE300temp
	    EndfE300=numpy.append(EndfE300,EndfE300temp)
	    barnsEndF300=numpy.append(barnsEndF300,barnsEndF300temp)
	    invEndfE300=numpy.append(invEndfE300,invEndfE300temp)
	
    log_printf(INFO,'Read in ENDF/B-VII XS correctly')

    #Plot values on top of each other
    fig=plt.figure()
    plt.loglog(EnT,barnsT)
    plt.loglog(EndfE300,barnsEndF300)
    CXStype=type_xs.title()
    plt.legend(['Doppler broadened '+El+'-'+A+' '+CXStype+' XS at temp=' + \
			temp,xssource+' '+El+'-'+A+' '+CXStype+' XS at temp=300K'], \
                                        loc='lower left',prop={'size':10})
    plt.grid()
    plt.title(CXStype+' Cross Section Comparison')
    plt.xlabel('XS [barns]')
    plt.ylabel('E [eV]')
    plt.savefig(dir+'/'+CXStype+'_XS_Comparison.png')
