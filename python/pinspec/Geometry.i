%module pinspec

%{
#include "../../src/Geometry.h"
#include "../../src/Region.h"
#define SWIG_FILE_WITH_INT
#include "../../src/Isotope.h"
#include "../../src/Material.h"
#include "../../src/Tally.h"
#include "../../src/Neutron.h"
#include "../../src/Fissioner.h"
#include "../../src/log.h"
#include "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
%}


%include ../../src/Geometry.h
%include ../../src/Region.h
%include ../../src/Isotope.h
%include ../../src/Material.h
%include ../../src/Tally.h
%include ../../src/Neutron.h
%include ../../src/Fissioner.h
%include ../../src/log.h
%include "numpy.i"

%init %{
    import_array();
%}

%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* xs, int n)}







