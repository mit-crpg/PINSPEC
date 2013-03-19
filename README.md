PINSPEC
=======

A Monte Carlo code for simple spectral calculations in nuclear reactor applications. To download the code, open
a terminal window and 'cd' into the directory where you would like to install PINSPEC. From here, download the code
as follows:

   > git clone https://github.com/wbinventor/PINSPEC.git PINSPEC
   
Now that you have downloaded the PINSPEC source code, it is time to install it. Enter the PINSPE directory:

   > cd PINSPEC

For a standard build and installation as a standalone Python package accessible to all users of your machine:

   > sudo python setup.py install
   
If you would rather install this as a standalone Python package only accessible to your username:

   > python setup.py install --user

Alternatively, you can do:

   > python2 setup.py install --prefix=$HOME/.local

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
