#! /usr/bin/env python

# System imports
import os
from distutils.sysconfig import get_config_vars
from distutils.core import *
from distutils      import sysconfig, dir_util

# Third-party modules - we depend on numpy for everything
import numpy

data_dir = os.path.join(sys.prefix, "local/lib/python2.7/dist-package/pinspec")


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

(opt,) = get_config_vars('OPT')
os.environ['OPT'] = ' '.join(
    flag for flag in opt.split() if flag != '-Wstrict-prototypes'
)

# range extension module
pinspec = Extension('_pinspec',
                   include_dirs=[numpy_include],
                   sources=['pinspec/Geometry.i',
                            'src/log.cpp', 'src/xsreader.cpp', 
                            'src/Isotope.cpp', 'src/Material.cpp', 
                            'src/Neutron.cpp', 'src/Tally.cpp', 
                            'src/Fissioner.cpp', 'src/Region.cpp', 
                            'src/Timer.cpp', 'src/Surface.cpp', 
                            'src/Geometry.cpp'],
                   extra_compile_args=['-O3', '-fopenmp',
                                    '-march=native', '-ffast-math'],
                   extra_link_args=['-lstdc++', '-fopenmp', '-lgomp', '-fPI'],
                   language='c++',
                   swig_opts=['-c++'],
                   )

# NumyTypemapTests setup
setup(  name        = 'PINSPEC',
        description = 'A monte carlo code for pin cell spectral calculations in nuclear reactor applications',
        author      = 'Will Boyd',
        author_email = 'wboyd@mit.edu',
        url = 'https://github.com/wbinventor/PINSPEC',
        version     = '0.1',
        ext_modules = [pinspec],
        packages = ['pinspec'],
		package_data = {'pinspec': ['xs-lib/*.txt', 'xs-lib/BackupXS/*.txt']},
        )

