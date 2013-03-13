#Input_Jessica.py
import matplotlib.pyplot as plt
import os
import sys
import numpy

class Nuclide(object):
    def __init__(self, name, ao):
       self.name=name
       self.ao=float(ao)
       self.AbsXS=numpy.array([])
       self.AbsXSE=numpy.array([])
       self.FisXS=numpy.array([])
       self.FisXSE=numpy.array([])
       self.ScatXS=numpy.array([])
       self.ScatXSE=numpy.array([])
       self.A=None
       self.El=None
       self.getAbsXS() 
       self.getFisXS()
       self.getScatXS()       
    def getAnEl(self):
       self.El, self.A=self.name.split("-", 1)
       self.A=float(self.A)
    def getAbsXS(self):
       cur_dir=os.getcwd()
       abxs_file=cur_dir+"/../XS_Lib/"+self.name+"-capture.txt"
       if os.path.exists(abxs_file)==True:
          print "Found "+self.name+" absorption cross section file"
       else:
          print "Did not find "+self.name+" absorption cross section file"
    #Parse abxs_file to read in line by line, skipping the first, etc.
       with open(abxs_file) as abxs:
          next(abxs)
          for line in abxs:
             E, XS=line.split('  ', 1)
             E=float(E)
             XS=float(XS)
             self.AbsXS=numpy.append(self.AbsXS,XS)
             self.AbsXSE=numpy.append(self.AbsXSE,E)
    def getFisXS(self):
       cur_dir=os.getcwd()
       fisxs_file=cur_dir+"/../XS_Lib/"+self.name+"-fission.txt"
       if os.path.exists(fisxs_file)==True:
          print "Found "+self.name+" fission cross section file"
       else:
          print "Did not find "+self.name+" fission cross section file"
    #Parse fisxs_file to read in line by line, skipping the first, etc.
       with open(fisxs_file) as fisxs:
          next(fisxs)
          for line in fisxs:
             E, XS=line.split('  ', 1)
             E=float(E)
             XS=float(XS)
             self.FisXS=numpy.append(self.FisXS,XS)
             self.FisXSE=numpy.append(self.FisXSE,E)
    def getScatXS(self):
       cur_dir=os.getcwd()
       scatxs_file=cur_dir+"/../XS_Lib/"+self.name+"-elastic.txt"
       if os.path.exists(scatxs_file)==True:
          print "Found "+self.name+" scattering cross section file"
       else:
          print "Did not find "+self.name+" scattering cross section file"
    #Parse scatxs_file to read in line by line, skipping the first, etc.
       with open(scatxs_file) as scatxs:
          next(scatxs)
          for line in scatxs:
             E, XS=line.split('  ', 1)
             E=float(E)
             XS=float(XS)
             self.ScatXS=numpy.append(self.ScatXS,XS)
             self.ScatXSE=numpy.append(self.ScatXSE,E)
    def plotXS(self):
       plt.loglog(self.AbsXSE,self.AbsXS)
       plt.loglog(self.FisXSE,self.FisXS)
       plt.loglog(self.ScatXSE,self.ScatXS)
       plt.legend(["Absorption","Fission", "Scattering"])
       plt.grid()
       plt.title('Cross section')
       plt.xlabel('XS [barns]')
       plt.ylabel('E [eV]')
       plt.show()
       plt.savefig("Cross_Sections.png")
