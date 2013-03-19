#Resonance Integral
import matplotlib.pyplot as plt
import os
import numpy
import math
import pinspec
#from scipy.integrate import simps

def compareXS(nameoffile, XStype='capture')
	#Resonance Integral Boundaries
	RIb=numpy.array([[0.01,0.1],[0.1,1.0],[6,10],[1,6],[10,25],[25,50],[50,100],[0.5,10000]], dtype=float)
	
	#Get fake XS from info given
	El, A, Rest = nameoffile.split('-', 2)
	#Find proper filename for fake XS
	if XStype=='scatter':
		XStype='elastic'
	path=str(pinspec.getXSLibDirectory())+'/'+El+'-'+A+'-'+XStype+'.txt')
	#Read in array for fictitious XS at 300
	#resT=open(path, 'r').readlines()
	EnT=numpy.array([])
	barnsT=numpy.array([])
	invEnT=numpy.array([])
	with open(path) as resT:
		#Parse out String containing Temperature
		Junk, T = resT.readline().split('=', 1)
		for line in resT:
			EnTt, barnsTt=line.split(',', 1)
			EnTt=float(EnTt)
			barnsTt=float(barnsTt)
			invEnTt=1/EnTt
			EnT=numpy.append(EnT,EnTt)
			barnsT=numpy.append(barnsT,barnsTt)
			invEnT=numpy.append(invEnT, invEnTt)
	print "Read in Doppler Broadened XS correctly"

	#Read in array for ENDF7 XS at 300
	npath=str(pinspec.getXSLibDirectory())+'/BackupXS/'+El+'-'+A+'-'+XStype+'.txt')
	#Endf300=open(npath, 'r').readlines()
	EndfE300=numpy.array([])
	barnsEndF300=numpy.array([])
	invEndfE300=numpy.array([])
	with open(npath) as Endf300:
		#Parse out String containing Temperature
		Junk, xssource = Endf300.readline().split(' ', 1)
		for line in Endf300:
			EndfE300t, barnsEndF300t=line.split(',', 1)
			EndfE300t=float(EndfE300t)
			barnsEndF300t=float(barnsEndF300t)
			invEndfE300t=1/EndfE300t
			EndfE300=numpy.append(EndfE300,EndfE300t)
			barnsEndF300=numpy.append(barnsEndF300,barnsEndF300t)
			invEndfE300=numpy.append(invEndfE300,invEndfE300t)
	print "read in ENDF/B-VII XS ok"

	#Plot values on top of each other
	fig=plt.figure()
	plt.loglog(EnT,barnsT)
	plt.loglog(EndfE300,barnsEndF300)
	plt.legend(["Doppler broadened "+El+'-'+A+' '+XStype+" XS at T="+T,xssource+' '+El+'-'+A+" at T=300K"])
	plt.grid()
	plt.title('Cross Section Comparison')
	plt.xlabel('XS [barns]')
	plt.ylabel('E [eV]')
	plt.savefig("XS_Comparison.png")

	#Make loop for RI calculation for Fake U238
	prod=numpy.zeros_like(barnsT)
	prodr=numpy.zeros_like(barnsEndF300)
	RIfake=numpy.zeros(len(RIb), dtype=float)
	RIreal=numpy.zeros(len(RIb), dtype=float)

	for i in range(len(RIb-1)):
		#Retrieve bounds
		Elow=RIb[i,0]
		Eupp=RIb[i,1]
		
		#Find index matching boundary in Energy vectors
		indlow=numpy.flatnonzero(EnT>=Elow) 
		indlo=numpy.array([indlow[0]], dtype=int)
		indupp=numpy.flatnonzero(EnT>=Eupp)
		indup=numpy.array([indupp[0]], dtype=int)
		indlowr=numpy.flatnonzero(EndfE300>=Elow) 
		indlor=numpy.array([indlowr[0]], dtype=int)
		induppr=numpy.flatnonzero(EndfE300>=Eupp)
		indupr=numpy.array([induppr[0]], dtype=int)

		#create vector to integrate, in	tegrate
		prod=(barnsT*invEnT)
		prodr=(barnsEndF300*invEndfE300)
		#en vector=EnT
		RIfake[i]=numpy.trapz(prod[indlo:indup],EnT[indlo:indup])
		RIreal[i]=numpy.trapz(prodr[indlor:indupr],EndfE300[indlor:indupr])
	print "Resonance Integral Bounds:"
	print RIb
	print "Doppler broadened "+El+'-'+A+' '+XStype+" Resonance Integrals at T="+T
	print str(RIfake)
	print "Real Resonance Integrals integrated from ENDF7 XS at T=300K"
	print str(RIreal)
		
