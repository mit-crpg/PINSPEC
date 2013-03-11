#Test array concatenation and printing
import numpy
import os

E=numpy.array([1,2,3,4,5,6])
XS=numpy.array([5,4,3,2,1,1])

EXS=(E,XS)
EXS=numpy.transpose(EXS)
cur_dir=os.getcwd()
out_name=cur_dir+"/../XS_Lib/U-238-capturef.txt"
numpy.savetxt(out_name,EXS,newline='\n',delimiter='  ')
#add Header
f=open(out_name)
text = f.read()
f.close()
# open the file again for writing
f = open(out_name, 'w')
f.write("U238 Doppler Broadened SLBW fictitious XS\n")
# write the original contents
f.write(text)
f.close()

#header='U238 Doppler Broadened SLBW fictitious XS'
