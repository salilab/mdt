%module _mdt
%include "typemaps.i"

%{
#include "../src/mdt.h"
#include "error_types.h"
#include "mod_error.h"
%}

%include "mdt_exceptions.i"
%include "mdt_lists.i"

typedef int mbool;

%apply (const float VARLIST[], int N_VARLIST) { (const float dx_dy[], int n_dx_dy) };

%include "../src/mdt.h"

%init {
#ifdef SWIGPYTHON
  mdterror = PyErr_NewException("_mdt.error", NULL, NULL);
  Py_INCREF(mdterror);
  PyModule_AddObject(m, "error", mdterror);
#endif
}
