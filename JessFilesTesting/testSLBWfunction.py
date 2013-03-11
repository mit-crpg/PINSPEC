#script that calls SLBWfunction

from SLBW import *

#user input
nameoffile = 'U-238-ResonanceParameters.txt'  # Must be Reich-Moore parameters
T=300 #Temp in Kelvin of target nucleus

#call function
SLBWXS(nameoffile,T)

