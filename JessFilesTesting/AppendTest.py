import os

class Nuclide(object):
    def __init__(self, name):
       self.name=name
       self.AbsXS=[]
       self.A=[]
       self.El=[]
    def getAnEl(self):
       self.El, self.A=self.name.split("-", 1)
       #Check whether this should be float or double
       self.A=float(self.A)
    def getAbsXS(self):
       cur_dir=os.getcwd()
       abxs_file=cur_dir+"/XS_Lib/"+self.name+"-capture.txt"
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

             
