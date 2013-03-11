#Resonance Integral
import matplotlib.pyplot as plt
import os
import numpy
import math
#from scipy.integrate import simps

#Read in array for fictitious XS at 300
cur_dir = os.getcwd()+'/U-238-capture300.txt'
res300=open(cur_dir, 'r').readlines()
En300=numpy.array([])
barns300=numpy.array([])
invEn300=numpy.array([])
for line in res300:
	En300t, barns300t=line.split('  ', 1)
        En300t=float(En300t)
	barns300t=float(barns300t)
	invEn300t=1/En300t
        En300=numpy.append(En300,En300t)
        barns300=numpy.append(barns300,barns300t)
	invEn300=numpy.append(invEn300, invEn300t)
print "read in fake XS ok"

#Read in array for ENDF7 XS at 300
cur_p = os.getcwd()+'/U-238-captureRealEndF7.txt'
Endf300=open(cur_p, 'r').readlines()
EndfE300=numpy.array([])
barnsEndF300=numpy.array([])
invEndfE300=numpy.array([])
for line in Endf300:
	EndfE300t, barnsEndF300t=line.split(',', 1)
	EndfE300t=float(EndfE300t)
	barnsEndF300t=float(barnsEndF300t)
        invEndfE300t=1/EndfE300t
        EndfE300=numpy.append(EndfE300,EndfE300t)
	barnsEndF300=numpy.append(barnsEndF300,barnsEndF300t)
	invEndfE300=numpy.append(invEndfE300,invEndfE300t)
print "read in ENDF7 XS ok"

#Plot values on top of each other
fig=plt.figure()
plt.loglog(En300,barns300)
plt.loglog(EndfE300,barnsEndF300)
plt.legend(["Fake U238 XS, 300","ENDFBVII U238 XS, 300"])
plt.grid()
plt.title('U238 Cross Section Comparison')
plt.xlabel('XS [barns]')
plt.ylabel('E [eV]')
plt.savefig("XS_Comparison_300.png")

#Resonance Integral Boundaries, can be passed in
RIb=numpy.array([[1e-5,1],[6,10],[1,6],[10,25],[25,50],[50,100],[100,1000]], dtype=float)

#Make loop for RI calculation for Fake U238
prod=numpy.zeros_like(barns300)
RI=numpy.zeros(len(RIb), dtype=float)

# Create an numpy interpolator object for the xs
#RI_interp_xs_object = interp1D(En300, barns300)
leg=[]
fig=plt.figure()
for i in range(len(RIb-1)):

	#Retrieve bounds
	Elow=RIb[i,0]
	Eupp=RIb[i,1]
	
	# Create array of energies between bounds
	#num_energies = 10
	#RI_energies = numpy.linspace(Elow, Eupp, num_energies)
	#RI_interp_xs = numpy.interp(RI_energies, En300, barns300)
	#plt.loglog(RI_energies,RI_interp_xs)
	#leg.append(str(i))
	#plt.loglog(En300, barns300)
	#RI[i]=numpy.trapz(RI_interp_xs, RI_energies)
	#if i==1:
		#print str(RI_energies)
		#print str(RI_interp_xs)
		#print str(Elow)
		#print str(Eupp)
	#Find index matching boundary in Energy vectors
	indlow=numpy.flatnonzero(En300>=Elow) 
	indlo=numpy.array([indlow[0]], dtype=int)
	indupp=numpy.flatnonzero(En300>=Eupp)
	indup=numpy.array([indupp[0]], dtype=int)
	#create vector to integrate, in	tegrate
	prod=(barns300*invEn300)
	#en vector=En300
	RI[i]=numpy.trapz(prod[indlo:indup],En300[indlo:indup])
	#RI[i]=prod[indlo:indup].sum()
	#RIt=binprod[:].sum()
print str(RI)
plt.legend(leg)
#plt.show()
	
