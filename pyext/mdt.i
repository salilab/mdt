%module _mdt
%include "typemaps.i"

%{
#include <glib.h>
#include "../src/mdt.h"
#include "../src/mdt_feature.h"
#include "../src/mdt_all_features.h"
#include "modeller.h"
%}

/* Ignore shared object import/export stuff */
#define MDTDLLEXPORT
#define MDTDLLLOCAL

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
%include "../src/mdt_feature.h"
%include "../src/mdt_all_features.h"

# Wrap MDT types
%include "mdt_type.i"
%include "mdt_library.i"

%init {
#ifdef SWIGPYTHON
  mdterror = PyErr_NewException("_mdt.MDTError", NULL, NULL);
  Py_INCREF(mdterror);
  PyModule_AddObject(m, "MDTError", mdterror);

  file_format_error = PyErr_NewException("_mdt.FileFormatError", mdterror,
                                         NULL);
  Py_INCREF(file_format_error);
  PyModule_AddObject(m, "FileFormatError", file_format_error);
#endif
}
