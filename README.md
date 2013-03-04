PINSPEC
=======

A Monte Carlo code for simple spectral calculations in nuclear reactor applications.

For a standard build, do: 
1) Comment out line 74 in src/Isotope.cpp; run 
   > swig -python -c++ SWIG/Region.i
   Then uncomment line 74. 

2) Run the build script,
   > sh build.sh
   Advanced options: --release, --debug, --profile, --benchmark 

3) cd into SWIG, do 
   > python
   > from region import *
   > fuel = Region()
   > ...

