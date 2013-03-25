import _pinspec
from pinspec import *
import os
import random

# Set a default logging level
log_setlevel(NORMAL)

# Set a log file name using a random integer id
setLogfileName('log/pinspec-' + str(random.randint(1,10000)) + '.log');

# Set the path to the xs library to the one that was installed
pkg_path = os.path.dirname(__file__)
xs_lib_path = os.path.join(pkg_path, 'xs-lib/')
setXSLibDirectory(xs_lib_path)

# Restore the cross-section library from backup
restoreXSLibrary()

# Create instances of TallyFactory and TallyBank singleton classes
pinspec.TallyFactory = TallyFactory.Get()
pinspec.TallyBank = TallyBank.Get()

