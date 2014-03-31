#! /usr/bin/env python

# System imports
import os
from distutils.sysconfig import get_config_vars
from distutils.core import *
from distutils.command.build_py import build_py
from distutils      import sysconfig, dir_util

# Third-party modules - we depend on numpy for everything
import numpy


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

(opt,) = get_config_vars('OPT')
os.environ['OPT'] = ' '.join(
    flag for flag in opt.split() if flag != '-Wstrict-prototypes'
)

#, '-L/opt/local/lib/gcc47'
#, '/opt/local/lib/gcc47/gcc/x86_64-apple-darwin11/4.7.2/include'

# range extension module
pinspec = Extension('_pinspec',
                   include_dirs=[numpy_include],
                   sources=['pinspec/Geometry.i',
                            'pinspec/src/log.cpp', 
                            'pinspec/src/xsreader.cpp', 
                            'pinspec/src/Isotope.cpp', 
                            'pinspec/src/Material.cpp', 
                            'pinspec/src/Neutron.cpp',
                            'pinspec/src/Tally.cpp',
                            'pinspec/src/TallyFactory.cpp',
			    'pinspec/src/TallyBank.cpp',
                            'pinspec/src/Fissioner.cpp',
                            'pinspec/src/Region.cpp',
                            'pinspec/src/Timer.cpp',
                            'pinspec/src/Surface.cpp', 
                            'pinspec/src/Geometry.cpp'],
                   extra_compile_args=['-O3', '-fopenmp',
                                    '-march=native', '-ffast-math', '-g'],
                   extra_link_args=['-lstdc++', '-fopenmp', '-lgomp'],
                   language='c++',
                   swig_opts=['-c++'],
                   )

# NumyTypemapTests setup
dist = setup(  name        = 'PINSPEC',
        description = 'A monte carlo code for pin cell spectral calculations in nuclear reactor applications',
        author      = 'Will Boyd',
        author_email = 'wboyd@mit.edu',
        url = 'https://github.com/wbinventor/PINSPEC',
        version     = '0.1',
        ext_modules = [pinspec],
        packages = ['pinspec'],
		package_data = {'pinspec': ['xs-lib/*.txt', 'xs-lib/BackupXS/*.txt']},
        )

# Rerun the build_py to setup links for C++ extension modules created by SWIG
# This prevents us from having to install twice
build_py = build_py(dist)
build_py.ensure_finalized()
build_py.run()
