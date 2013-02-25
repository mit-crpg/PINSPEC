%module material
%{
#include "../src/Isotope.h"
#include "../src/Material.h"
%}

// Very simple C++ example to create materials and isotopes using SWIG

%include ../src/Isotope.h
%include ../src/Material.h

