#Parse string
import numpy
import os
import sys

#function to separate a value into number and exponent, and return the float
def convert(string):
    n=8
    num, exp=[string[i:i+n] for i in range(0, len(string), n)]
    num=float(num)
    exp=float(exp)
    value=num*(10**exp)
    return value

#function to parse a given string line into the individual float values of E0, GN, and GG for future appending
def parse(string):
    result=numpy.array([])
    E0, J, GN, GG, GFA, GFB=string.strip().split(' ', 5)
    #Further parse individual strings, convert to floats
    E0=convert(E0)
    GN=convert(GN)
    GG=convert(GG)
    together=E0, GN, GG
    result=numpy.append(result,together)
    return result

string1=' 6.673491+0 5.000000-1 1.475792-3 2.300000-2 0.000000+0 9.990000-9'
string2=' 2.087152+1 5.000000-1 1.009376-2 2.286379-2 5.420000-8 0.000000+0'

out1=parse(string1)
out2=parse(string2)
print out1
print out2

E0=numpy.array([])
E0=numpy.append(E0,out1[0])
E0=numpy.append(E0,out2[0])

GN=numpy.array([])
GN=numpy.append(GN,out1[1])
GN=numpy.append(GN,out2[1])

GG=numpy.array([])
GG=numpy.append(GG,out1[2])
GG=numpy.append(GG,out2[2])

print E0
print GN
print GG






