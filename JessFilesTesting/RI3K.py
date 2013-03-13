#Resonance Integral
import matplotlib.pyplot as plt
import os
import numpy
import math
#from scipy.integrate import simps

#Read in array for fictitious XS at 3
cur_dir = os.getcwd()+'/U-238-capture3.txt'
res3=open(cur_dir, 'r').readlines()
En3=numpy.array([])
barns3=numpy.array([])
invEn3=numpy.array([])
for line in res3:
	En3t, barns3t=line.split('  ', 1)
        En3t=float(En3t)
	barns3t=float(barns3t)
	invEn3t=1/En3t
        En3=numpy.append(En3,En3t)
        barns3=numpy.append(barns3,barns3t)
	invEn3=numpy.append(invEn3, invEn3t)
print "read in XS at T=3K ok"

#Read in array for our XS at 300
cur_p = os.getcwd()+'/U-238-capture300.txt'
res300=open(cur_p, 'r').readlines()
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
	invEn300=numpy.append(invEn300,invEn300t)
print "read in XS at T=300 ok"

#Plot values on top of each other
fig=plt.figure()
plt.loglog(En3,barns3)
plt.loglog(En300,barns300)
plt.legend(["Fake U238 XS at T=3K","Fake U238 XS at T=300K"])
plt.grid()
plt.title('U238 Cross Section Comparison at Different Temperatures')
plt.xlabel('XS [barns]')
plt.ylabel('E [eV]')
plt.savefig("XS_Comparison_3v300.png")

#Resonance Integral Boundaries, can be passed in
RIb=numpy.array([[0.01,0.1],[0.1,1.0],[6,10],[1,6],[10,25],[25,50],[50,100],[0.5,10000]], dtype=float)

#Make loop for RI calculation for Fake U238
prod=numpy.zeros_like(barns3)
prodr=numpy.zeros_like(barns300)
RIfake=numpy.zeros(len(RIb), dtype=float)
RIreal=numpy.zeros(len(RIb), dtype=float)

for i in range(len(RIb-1)):
	#Retrieve bounds
	Elow=RIb[i,0]
	Eupp=RIb[i,1]
	
	#Find index matching boundary in Energy vectors
	indlow=numpy.flatnonzero(En3>=Elow) 
	indlo=numpy.array([indlow[0]], dtype=int)
	indupp=numpy.flatnonzero(En3>=Eupp)
	indup=numpy.array([indupp[0]], dtype=int)
	indlowr=numpy.flatnonzero(En300>=Elow) 
	indlor=numpy.array([indlowr[0]], dtype=int)
	induppr=numpy.flatnonzero(En300>=Eupp)
	indupr=numpy.array([induppr[0]], dtype=int)

	#create vector to integrate, in	tegrate
	prod=(barns3*invEn3)
	prodr=(barns300*invEn300)
	#en vector=En3
	RIfake[i]=numpy.trapz(prod[indlo:indup],En3[indlo:indup])
	RIreal[i]=numpy.trapz(prodr[indlor:indupr],En300[indlor:indupr])
print "Resonance Integral Bounds:"
print RIb
print "Fake U238 Resonance Integrals at T=3"
print str(RIfake)
print "Fake U238 Resonance Integrals at T=300"
print str(RIreal)
	
