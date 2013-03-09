%module pinspec

%{
    #define SWIG_FILE_WITH_INIT
    #include "../../src/Geometry.h"
    #include "../../src/Region.h"
    #include "../../src/Isotope.h"
    #include "../../src/Material.h"
    #include "../../src/Tally.h"
    #include "../../src/Neutron.h"
    #include "../../src/Fissioner.h"
    #include "../../src/log.h"
%}

%include "numpy.i"

%init %{
     import_array();
%}

%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* xs, int num_xs)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* energies, int num_xs)}

%include ../../src/Geometry.h
%include ../../src/Region.h
%include ../../src/Isotope.h
%include ../../src/Material.h
%include ../../src/Tally.h
%include ../../src/Neutron.h
%include ../../src/Fissioner.h
%include ../../src/log.h
