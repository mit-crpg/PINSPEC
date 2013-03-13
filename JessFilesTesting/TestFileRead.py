#Input file read methods
import os
import numpy

#user input
nameoffile='Pu-241-ENDF.txt' #Must be Reich-Moore parameters
numberofpositiveresonances=14

#Get file path
cur_dir=os.getcwd()
filepath=cur_dir+"/../XS_Lib/"+nameoffile
if os.path.exists(filepath)==True:
   print "Found resonance parameter file"
else:
   print "Did not find resonance parameter file"
#Initialize desired arrays
E0=numpy.array([])
GN=numpy.array([])
GG=numpy.array([])
restxt=open(filepath, 'r').readlines()

for line in restxt:
   if '4*pi' in line:
      #Parse first line for SigP
      junk, SigP, barns=line.split(' ', 2)
      SigP=float(SigP)
      break
for line in restxt:
   if '  ----------' in line:
      gib=line
      print 'hey' 
##try a counter or something. This sucks. 

switch=False
for line in restxt:
   if switch==True:
      #See if first line is negative
      num, restofline =line.split('.', 1)
      num=float(num)
      #skip negative resonances
      if num<0:
         next(restxt)
      else:
         print line
      if len(data)>1:
          print line
      else:
          break
   if '  ----------' in line:
      switch=True
    

