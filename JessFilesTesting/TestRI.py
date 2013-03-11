#Resonance Integral
import matplotlib.pyplot as plt
import os
import numpy
import math


#Read in array for ENDF7 XS at 300
cur_p = os.getcwd()+'/U-238-captureRealEndF7.txt'
Endf300=open(cur_p, 'r').readlines()
EndfE300=numpy.empty_like(Endf300)
barnsEndF300=numpy.empty_like(Endf300)
for line in Endf300:
	EndfE300t, barnsEndF300t=line.split(',', 1)
	EndfE300t=float(EndfE300t)
	barnsEndF300t=float(barnsEndF300t)
        EndfE300=numpy.append(EndfE300,EndfE300t)
	barnsEndF300=numpy.append(barnsEndF300,barnsEndF300t)
