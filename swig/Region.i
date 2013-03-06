%module region
%{
#include "../src/Region.h"
#include "../src/Isotope.h"
#include "../src/Material.h"
#include "../src/Tally.h"
#include "../src/Neutron.h"
%}

// Very simple C++ example to create regions using SWIG

%include ../src/Region.h
%include ../src/Isotope.h
%include ../src/Material.h
%include ../src/Tally.h
%include ../src/Neutron.h

