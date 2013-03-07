%module pinspec

%{
#include "../../src/Geometry.h"
#include "../../src/Region.h"
#include "../../src/Isotope.h"
#include "../../src/Material.h"
#include "../../src/Tally.h"
#include "../../src/Neutron.h"
#include "../../src/Fissioner.h"
#include "../../src/log.h"
%}


%include ../../src/Geometry.h
%include ../../src/Region.h
%include ../../src/Isotope.h
%include ../../src/Material.h
%include ../../src/Tally.h
%include ../../src/Neutron.h
%include ../../src/Fissioner.h
%include ../../src/log.h



%{
	#define SWIG_FILE_WITH_INT
	#include "../../src/Isotope.h"
	
%}

%include "numpy.i"

%init %{
	import_array();
%}

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* xs, int n)}

%include "../../src/Isotope.h" 










