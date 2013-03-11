#SLBW.py
import matplotlib.pyplot as plt
import os
import sys
import numpy
import math
import scipy.special as spec #newest version of scipy,needed Cython 0.18 was needed, installed using python-pip, command pip install cython, pulled scipy from github

#User defined parameters
nameoffile='U-238-ENDF.txt' #Must be Reich-Moore parameters
numberofpositiveresonances=14
Emin=1e-5            #Min energy of fictitious SLBW XS
Ecut=1000            #Max energy of fictitious SLBW XS
Ebnwdth=0.075        #Resolution on the energy scale for the cross sections
T=293.15             #Temp in Kelvin of target nucleus
resspacing=25        #resonance spacing in eV for fake identical resonances
idntclfakereslb=300  #Lower bound of fake identical resonances after SLBW; THIS SHOULD BE AUTOMATED
idntclfakeresub=1000 #Upper cutoff of identical fake resonances before unresolved region
gamg=0.023           #Gamma Gamma values for identical fake resonances

#---------------------------------------------

#Get file path
cur_dir=os.getcwd()
filepath=cur_dir+"/../XS_Lib/"+nameoffile
if os.path.exists(filepath)==True:
   print "Found resonance parameter file"
else:
   print "Did not find resonance parameter file"

restxt=open(filepath, 'r').readlines()
for line in restxt:
   if '4*pi' in line:
      #Parse first line for SigP
      junk, SigP, barns=line.split(' ', 2)
      SigP=float(SigP)
      break
x=1
for line in restxt:
      if '  ----------' in line:
         print x
         break
      else:
         x+=1
y=0
while y<x:
   next(restxt) #not in file anymore
   y+=1
num, restofline =restxt.readline.split('.', 1)
num=float(num)
print num

