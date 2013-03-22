from pinspec import *
import _pinspec
import os

# Set the path to the xs library to the one that was installed
pkg_path = os.path.dirname(__file__)
xs_lib_path = os.path.join(pkg_path, 'xs-lib/')
setXSLibDirectory(xs_lib_path)

pinspec.TallyFactory = TallyFactory.Get()
pinspec.TallyBank = TallyBank.Get()
