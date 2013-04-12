PINSPEC
=======

A Monte Carlo code for simple spectral calculations in nuclear reactor applications. To download the code, open
a terminal window and 'cd' into the directory where you would like to install PINSPEC. From here, download the code
as follows:

   > git clone https://github.com/wbinventor/PINSPEC.git PINSPEC
   
Now that you have downloaded the PINSPEC source code, it is time to install it. Enter the PINSPEC directory:

   > cd PINSPEC

For a standard build and installation as a standalone Python package accessible from any directory on your machine:

   > python setup.py install --user

To access the C++ classes from a Python interpreter, do:

   > python

   > from pinspec import *

   > fuel = Region()

   > ...
