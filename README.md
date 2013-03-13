PINSPEC
=======

A Monte Carlo code for simple spectral calculations in nuclear reactor applications.

For a standard build and installation as a standalone Python package:

   > cd PINSPEC

   > sudo python setup.py install

To run a sample infinite medium input file, do: 

   > cd sample-input

   > python infinite.py

To run a sample homogeneouse equivalence input file, do: 

   > cd sample-input

   > python equivalence.py

To access the C++ classes in Python, do

   > python

   > from pinspec import *

   > fuel = Region()

   > ...
