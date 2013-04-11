import unit_tests
import sample_geometries
import sample_problems
import unittest
from pinspec import *
from pinspec.log import *

##########################################
##########  UNIT TESTING SUITE ###########
##########################################
# 1) Unit test important C++ and python functions
# 2) Run sample geometries and generate sample output
# 3) Run test problems to showcase capabilities to code

# set debug level to UNITTEST
py_setlevel('UNITTEST')

# create unittest loader and test suite
loader = unittest.TestLoader()
test_suite = unittest.TestSuite()

# 1) Unit test important C++ and python functions
test_suite.addTests(loader.loadTestsFromModule(unit_tests))

# 2) Run sample geometries and generate sample output
test_suite.addTests(loader.loadTestsFromModule(sample_geometries))

# 3) Run test problems to showcase capabilities to code
test_suite.addTests(loader.loadTestsFromModule(sample_problems))

# Run unit tests
unittest.TextTestRunner().run(test_suite)
