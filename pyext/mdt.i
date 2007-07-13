%module _mdt
%include "typemaps.i"

%{
#include <glib.h>
#include "../src/mdt.h"
#include "mod_error_types.h"
#include "modeller.h"
#include "mod_error.h"
%}

%include "fortran-pointers.i"
%include "mdt_exceptions.i"
%include "mdt_lists.i"
%include "mdt_glib.i"

%apply double *OUTPUT { double * };
%apply (const float VARLIST[], int N_VARLIST) { (const float dx_dy[], int n_dx_dy) };
%apply (const int VARLIST[], int N_VARLIST) { (const int features[], int n_features) };
%apply (const int VARLIST[], int N_VARLIST) { (const int offset[], int n_offset) };
%apply (const int VARLIST[], int N_VARLIST) { (const int shape[], int n_shape) };
%apply (const int VARLIST[], int N_VARLIST) { (const int indices[], int n_indices) };

%include "../src/mdt.h"
%include "../src/mdt_alignment.h"

# Wrap MDT types
%include "mdt_type.i"
%include "mdt_library.i"

%init {
#ifdef SWIGPYTHON
  mdterror = PyErr_NewException("_mdt.error", NULL, NULL);
  Py_INCREF(mdterror);
  PyModule_AddObject(m, "error", mdterror);
#endif
}
