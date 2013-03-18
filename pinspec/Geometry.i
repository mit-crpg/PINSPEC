%module pinspec

%{
    #define SWIG_FILE_WITH_INIT
    #include "../src/Geometry.h"
    #include "../src/Region.h"
    #include "../src/Isotope.h"
    #include "../src/Material.h"
    #include "../src/Tally.h"
    #include "../src/Neutron.h"
    #include "../src/Fissioner.h"
    #include "../src/log.h"
    #include "../src/xsreader.h"
%}

%include "numpy.i"

%init %{
     import_array();
%}

%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* xs, int num_xs)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* energies, int num_xs)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* data, int num_bins)}
%apply (double* IN_ARRAY1, int DIM1) {(double* edges, int num_edges)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* cdf, int num_bins)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* cdf_energies, int num_bins)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* cdfs, int num_values)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* Eprime_to_E, int num_bins)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* E_to_kT, int num_cdfs)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* dist, int num_values)}


%include ../src/Geometry.h
%include ../src/Region.h
%include ../src/Isotope.h
%include ../src/Material.h
%include ../src/Tally.h
%include ../src/Neutron.h
%include ../src/Fissioner.h
%include ../src/log.h
%include ../src/xsreader.h

