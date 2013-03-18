# SLBW.py
import pinspec
import matplotlib.pyplot as plt
import numpy
import os
import shutil
import math
import scipy.special as spec  # newest version of scipy,needed Cython 0.18 was needed, installed using python-pip, command pip install cython, pulled scipy from github

def SLBWXS(nameoffile,T,typeXS='capture',numberofpositiveresonances = 14,Emin = 1e-5,Ecut = 1000.0,Ebnwdth = 0.075,resspacing = 25,idntclfakereslb = 300,gamg = 0.023,flatxs = 0.1,finalEcut=20E6):
#---------------------------------------------------------------------------------
# Function that creates resonant capture and scatter XS at a given temperature and
# frankensteins it with given information
#----------------------------------------------------------------------------------
        # User defined parameters
	#nameoffile = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
	#numberofpositiveresonances = 14
	#Emin = 1e-5  # Min energy of fictitious SLBW XS
	#Ecut = 1000.0  # Max energy of fictitious SLBW XS
	#Ebnwdth = 0.075  # Resolution on the energy scale for the cross sections
	#T = 300  # Temp in Kelvin of target nucleus
	#resspacing = 25  # resonance spacing in eV for fake identical resonances
	#idntclfakereslb = 300  # Lower bound of fake identical resonances before unresolved (upper bound is Ecut)
	#gamg = 0.023  # Gamma Gamma values for identical fake resonances
	#flatxs = 0.1
	#finalEcut=20E6


	#---------------------------------------------
	#---Parsing Functions
	#---------------------------------------------
	# function to separate a value into number and exponent, and return the float
	def convert(string):
		n = 8
		num, exp = [string[i:i + n] for i in range(0, len(string), n)]
		num = float(num)
		exp = float(exp)
		value = num * (10 ** exp)
		return value
	#----------------------------------------------
	# function to parse a given res param string line into the individual float values of E0, GN, and GG for future appending
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

	#-------------------------------------------------
	#---Construct resonance parameter arrays from given file
	#--------------------------------------------------
	# Get file path
	#cur_dir = os.getcwd()
	filepath = str(pinspec.getXSLibDirectory()) + nameoffile
	if os.path.exists(filepath) == True:
		pinspec.log_printf(pinspec.INFO, 'Loading resonance paramater file' + filepath)
	else:
		pinspec.log_printf(pinspec.WARNING, 'Unable to load resonance parameter file' + filepath)
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
		# parse each of the rest of the desired resonances and append to the arrays
		for i, line in enumerate(restxt):
			if i < (numberofpositiveresonances - 1):
				out = parse(line)
				E0 = numpy.append(E0, out[0])
				GN = numpy.append(GN, out[1])
				GG = numpy.append(GG, out[2])
			else:
				break
	restxt.close()
	#--------------------------------------------------
	# Now we have SigP, E0, GN, and GG for first 14 (for U238) resonances, need  rest
	#--------------------------------------------------
	#resonances up to 1keV------------
	# Make up an E_0 vector with resonance spacing from 300 to 1 keV
	numbns = round((Ecut - idntclfakereslb) / resspacing)
	Enot = numpy.linspace(idntclfakereslb, Ecut, numbns)

	# Make up Gamma_gamma (0.023 for U-238)
	Gamma_g = numpy.empty_like(Enot)
	Gamma_g[:] = gamg

	# Make up gamma_n for each
	E_last = E0[numberofpositiveresonances - 1]  # Last known resonance E_0 (291.0206 for U-238)
	Gamma_nlast = GN[numberofpositiveresonances - 1]  # Last known Gn value
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

	# pull A value from filename
	El, A, rest = nameoffile.split('-', 2)
	A = float(A)

	# Construct Energy grid for fictitious XS
	nbins = round((Ecut - Emin) / Ebnwdth)
	E = numpy.logspace(numpy.log10(Emin), numpy.log10(Ecut), nbins)

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
	squig = Gamma * ((A / (4 * k * T * E0)) ** 0.5)  # squiggle
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
		#plt.plot(psi)
		#plt.show()
		chi = 2.0*psichi.imag
		#plt.plot(chi)
		#plt.show()

		# generate cross sections for this resonance, and sum
		sigma_g = sigma_g + (GN[i] / Gamma[i]) * (GG[i] / Gamma[i]) * ((E0[i] / E) ** 0.5) * (r[i] * psi)
		sigma_n = sigma_n + ((GN[i] / Gamma[i])**2) * (r[i] * psi + q[i] * chi)
		#plt.loglog(E, sigma_n)
		#plt.show()		
	
	#plot fictitious XS
	#plt.loglog(E, sigma_n)
	#plt.savefig('U238SXS.png')
	#-----------------------------------------------------------------
	# Now that we have cross sections, we can frankenstein after 1 keV
	#----------------------------------------------------------------
	# From 1keV to 20E6 eV, use 0.1 barns for capture (automated to use inputs)
	E2 = numpy.linspace(Ecut + Ebnwdth, finalEcut, 10)
	sigma_g2 = numpy.empty_like(E2)
	sigma_g2[:] = flatxs
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


	#plot fictitious XS
	#plt.loglog(E, sigma_n)
	#plt.savefig('U238SXS.png')
	#plot fictitious XS
	#plt.loglog(E, sigma_g)
	#plt.savefig('U238AXS.png')

	#----------------------------------
	#---So if capture is specified
	#---------------------------------
	if typeXS=='capture':
		# write output file for capture XS
		out_name = pinspec.getXSLibDirectory()+El+'-'+str(int(A))+'-capture.txt'
		numpy.savetxt(out_name, EAXS, newline='\n', delimiter=',')
		# go back and add header
		f = open(out_name)
		text = f.read()
		f.close()
		# open the file again for writing
		f = open(out_name, 'w')
		f.write('Doppler Broadened SLBW fictitious capture XS at T='+str(T)+'K\n')
		# write the original contents
		f.write(text)
		f.close()
	#----------------------------------
	#---So if capture is specified
	#---------------------------------
	if typeXS=='scatter':
		# write output file for scatter XS
		out_names = pinspec.getXSLibDirectory()+El+'-'+str(int(A))+'-elastic.txt'
		numpy.savetxt(out_names, ESXS, newline='\n', delimiter=',')
		# go back and add header
		g = open(out_names)
		texts = g.read()
		g.close()
		# open the file again for writing
		g = open(out_names, 'w')
		g.write('Doppler Broadened SLBW fictitious resonant scattering XS at T='+str(T)+'K\n')
		# write the original contents
		g.write(texts)
		g.close()
		#-----------------
def generatePotentialScattering(nameoffile,Emin = 1e-5,finalEcut=20E6):
#----------------------------------------------------------------------
# Function that creates  flat scattering XS at sigma_pot
#----------------------------------------------------------------------
	#nameoffile = 'U-238-ResonanceParameters.txt'
	filepath = str(pinspec.getXSLibDirectory()) + nameoffile
	if os.path.exists(filepath) == True:
		pinspec.log_printf(pinspec.INFO, 'Loading resonance paramater file' + filepath)
	else:
		pinspec.log_printf(pinspec.WARNING, 'Unable to load resonance parameter file' + filepath)
	with open(filepath) as restxt:
		# Parse first line for SigP
		junk, SigP, barns = restxt.readline().split(' ', 2)
		SigP = float(SigP)
	E = numpy.logspace(numpy.log10(Emin), numpy.log10(finalEcut), 2)
	SXS=numpy.zeros_like(E)
	SXS[:]=SigP
	ESXS = (E, SXS)
	ESXS = numpy.transpose(ESXS)
	# pull A value from filename
	El, A, rest = nameoffile.split('-', 2)
	A = float(A)
	# write output file for scatter XS
	out_names = pinspec.getXSLibDirectory()+El+'-'+str(int(A))+'-elastic.txt'
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
	#-----------------
#def replaceResXS():
#---------------------------------------------------------------------------
#Function that cleans up library by replacing XS files with ENDFBVII files
#---------------------------------------------------------------------------

