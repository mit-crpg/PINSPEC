PINSPEC
=======

A Monte Carlo code for simple spectral calculations in nuclear reactor applications.

For a standard build, do: 

   > sh build.sh

   of 

   > ./build.sh

   Advanced options: --release, --debug, --profile, --benchmark 

To use the infinite geometry class, do 

   > cd python

   > python2 infinite.py


To access the C++ classes in Python, do

   > cd swig

   > python

   > from region import *

   > fuel = Region()

   > ...

