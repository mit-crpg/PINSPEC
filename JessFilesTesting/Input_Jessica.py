#Input_Jessica.py
#import matplotlib.pyplot as plt
import os
import sys
import numpy

class Nuclide(object):
    def __init__(self, name, ao):
       self.name=name
       self.ao=float(ao)
       self.AbsXS=[]
       self.FisXS=[]
       self.ScatXS=[]
       self.A=[]
       self.El=[]
    def getAnEl(self):
       self.El, self.A=self.name.split("-", 1)
       #Check whether this should be float or double
       self.A=float(self.A)
    def getAbsXS(self):
       cur_dir=os.getcwd()
       abxs_file=cur_dir+"/../XS_Lib/"+self.name+"-capture.txt"
       if os.path.exists(abxs_file)==True:
          print "Found"+self.name+"absorption cross section file"
       else:
          print "Did not find"+self.name+" absorption cross section file"
    #Parse abxs_file to read in line by line, skipping the first, etc.
       with open(abxs_file) as abxs:
          next(abxs)
          for line in abxs:
             E, XS=line.split('  ', 1)
             E=float(E)
             XS=float(XS)
             #There must be a better way to do this and make it a real array
             data=E,XS
             self.AbsXS.append(data)
    def getFisXS(self):
       cur_dir=os.getcwd()
       fisxs_file=cur_dir+"/../XS_Lib/"+self.name+"-fission.txt"
       if os.path.exists(fisxs_file)==True:
          print "Found"+self.name+"fission cross section file"
       else:
          print "Did not find"+self.name+"fission cross section file"
    #Parse fisxs_file to read in line by line, skipping the first, etc.
       with open(fisxs_file) as fisxs:
          next(fisxs)
          for line in fisxs:
             E, XS=line.split('  ', 1)
             E=float(E)
             XS=float(XS)
             #There must be a better way to do this and make it a real array
             data=E,XS
             self.FisXS.append(data)
    def getScatXS(self):
       cur_dir=os.getcwd()
       
       scatxs_file=cur_dir+"/../XS_Lib/"+self.name+"-elastic.txt"
       if os.path.exists(scatxs_file)==True:
          print "Found"+self.name+"scattering cross section file"
       else:
          print "Did not find"+self.name+"scattering cross section file"
    #Parse scatxs_file to read in line by line, skipping the first, etc.
       with open(scatxs_file) as scatxs:
          next(scatxs)
          for line in scatxs:
             E, XS=line.split('  ', 1)
             E=float(E)
             XS=float(XS)
             #There must be a better way to do this and make it a real array
             data=E,XS
             self.ScatXS.append(data)
