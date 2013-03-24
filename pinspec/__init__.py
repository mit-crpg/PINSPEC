import _pinspec
from pinspec import *
import os

# Set the path to the xs library to the one that was installed
log_setlevel(NORMAL)

pkg_path = os.path.dirname(__file__)
xs_lib_path = os.path.join(pkg_path, 'xs-lib/')
setXSLibDirectory(xs_lib_path)

restoreXSLibrary()

pinspec.TallyFactory = TallyFactory.Get()
pinspec.TallyBank = TallyBank.Get()

