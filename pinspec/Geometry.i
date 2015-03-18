%module pinspec

%{
    #define SWIG_FILE_WITH_INIT
    #include "src/Geometry.h"
    #include "src/Region.h"
    #include "src/Surface.h"
    #include "src/Isotope.h"
    #include "src/Material.h"
    #include "src/Tally.h"
    #include "src/TallyBank.h"
    #include "src/TallyFactory.h"
    #include "src/Neutron.h"
    #include "src/Fissioner.h"
    #include "src/log.h"
    #include "src/vector.h"
    #include "src/xsreader.h"
    #include "src/Timer.h"

    #define printf PySys_WriteStdout

    /* Exception helpers */
    static int swig_c_error_num = 0;
    static char swig_c_err_msg[1024];

    const char* err_occurred(void) {
      if (swig_c_error_num) {
          swig_c_error_num = 0;
          return (const char*)swig_c_err_msg;
      }
      return NULL;
    }

    void set_err(const char *msg) {
      swig_c_error_num = 1;
      strncpy(swig_c_err_msg, msg, 1024);
    }
%}


%exception {
  try {
    $function
  } catch (const std::exception &e) {
    SWIG_exception(SWIG_RuntimeError, e.what()); 
  }
}


%include "numpy.i"


%init %{
     import_array();
%}


%ignore Tally::operator=(Tally* tally);
%ignore Tally::operator=(const Tally& tally);


%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* xs, int num_xs)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* energies, int num_xs)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* data, int num_bins)}
%apply (double* IN_ARRAY1, int DIM1) {(double* edges, int num_edges)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* cdf, int num_bins)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* cdf_energies, int num_bins)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* cdfs, int num_values)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* Eprime_to_E, int num_bins)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* E_to_kT, int num_cdfs)}
%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* pdfs, int num_values)}
%apply (double* IN_ARRAY1, int DIM1) {(double* energies, int num_energies), (double* elastic_xs, int num_xs)}
%apply (double* IN_ARRAY1, int DIM1) {(double* energies, int num_energies), (double* capture_xs, int num_xs)}
%apply (double* IN_ARRAY1, int DIM1) {(double* energies, int num_energies), (double* fission_xs, int num_xs)}


%apply (int* IN_ARRAY1, int DIM1) {(const int* amt, const int length)}
%apply (float* IN_ARRAY1, int DIM1) {(const float* amt, const int length)}
%apply (double* IN_ARRAY1, int DIM1) {(const double* amt, const int length)}

%include <exception.i>
%include src/Geometry.h
%include src/Region.h
%include src/Surface.h
%include src/Isotope.h
%include src/Material.h
%include src/Tally.h
%include src/TallyBank.h
%include src/TallyFactory.h
%include src/Neutron.h
%include src/Fissioner.h
%include src/log.h
%include src/vector.h
%include src/xsreader.h
%include src/Timer.h


#define printf PySys_WriteStdout
