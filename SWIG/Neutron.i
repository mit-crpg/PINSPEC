%module neutron
%{
#include "../src/Neutron.h"
%}

// Very simple C++ example with Neutron to create Neutron struct for SWIG

%include ../src/Neutron.h

%{
extern neutron* initializeNewNeutron();
%}

extern neutron* initializeNewNeutron();

