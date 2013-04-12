import _pinspec
from pinspec import *
import os
import random
import datetime

# Set a default logging level
log_setlevel(NORMAL)

# Set a log file name using a date and time
now = datetime.datetime.now()
current_time = str(now.month) + '-' + str(now.day) + '-' + str(now.year) + '--' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second)
setLogfileName('log/pinspec-' + current_time + '.log');

# Set the path to the xs library to the one that was installed
pkg_path = os.path.dirname(__file__)
xs_lib_path = os.path.join(pkg_path, 'xs-lib/')
setXSLibDirectory(xs_lib_path)

# Restore the cross-section library from backup
restoreXSLibrary()

# Get instances of TallyBank and TallyFactory singleton classes
pinspec.TallyBank = TallyBank.Get()
pinspec.TallyFactory = TallyFactory.Get()

# Get instances of SurfaceFactory and RegionFactory singleton classes
pinspec.SurfaceBank = SurfaceFactory.Get()
pinspec.RegionFactory = RegionFactory.Get()
