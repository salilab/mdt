%module mdt
%include "typemaps.i"

%{
#include <glib.h>
#include "mdt.h"
#include "mdt_feature.h"
#include "mdt_all_features.h"
#include "modeller.h"
%}

/* Ignore shared object import/export stuff */
#define MDTDLLEXPORT
#define MDTDLLLOCAL

/* Pass scripting language objects to new_foo() routines as opaque pointers */
%typemap(in) gpointer scriptobj {
  $1 = $input;
}

%include "fortran-pointers.i"
%include "mdt_exceptions.i"
%include "mdt_lists.i"
%include "mdt_glib.i"
%include "mdt_callbacks.i"

%apply double *OUTPUT { double * };
%apply (const float VARLIST[], int N_VARLIST) { (const float dx_dy[], int n_dx_dy) };
%apply (const int VARLIST[], int N_VARLIST) { (const int features[], int n_features) };
%apply (const int VARLIST[], int N_VARLIST) { (const int offset[], int n_offset) };
%apply (const int VARLIST[], int N_VARLIST) { (const int shape[], int n_shape) };
%apply (const int VARLIST[], int N_VARLIST) { (const int indices[], int n_indices) };
%apply (const int VARLIST[], int N_VARLIST) { (const int chunk_size[], int n_chunk_size) };

%include "mdt.h"
%include "mdt_alignment.h"
%include "mdt_all_features.h"

/* Wrap MDT types */
%include "mdt_type.i"
%include "mdt_library.i"
%include "mdt_feature.i"

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

%{
#ifdef MDT_WITH_NUMPY
#include <numpy/arrayobject.h>

static int import_numpy_module(void)
{
  static gboolean imported;

  if (imported) {
    return 0;
  } else {
    int ret = _import_array();
    if (ret == 0) {
      imported = TRUE;
    }
    return ret;
  }
}

static int get_numpy_bin_type(const struct mdt *mdt)
{
  switch (mdt->base.bin_type) {
  case MOD_MDTB_FLOAT:
    return NPY_FLOAT;
  case MOD_MDTB_DOUBLE:
    return NPY_DOUBLE;
  case MOD_MDTB_INT32:
    return NPY_INT32;
  case MOD_MDTB_UINT32:
    return NPY_UINT32;
  case MOD_MDTB_INT16:
    return NPY_INT16;
  case MOD_MDTB_UINT16:
    return NPY_UINT16;
  case MOD_MDTB_INT8:
    return NPY_INT8;
  case MOD_MDTB_UINT8:
    return NPY_UINT8;
  }
  g_assert_not_reached();
  return 0;
}
#endif
%}


%inline %{
PyObject *get_numpy(struct mdt *mdt, PyObject *mdt_pyobj)
{
#ifdef MDT_WITH_NUMPY
  int i, type_num;
  PyObject *obj;
  npy_intp *dims;
  char *data = mdt->base.bindata;

  if (import_numpy_module() != 0) {
    return NULL;
  }

  dims = g_malloc(sizeof(npy_intp) * mdt->base.nfeat);
  for (i = 0; i < mdt->base.nfeat; ++i) {
    dims[i] = mdt->base.features[i].nbins;
  }

  type_num = get_numpy_bin_type(mdt);
  /* Note that MDT tables are C-style contiguous so no special strides or
     other flags need to be passed to NumPy */
  obj = PyArray_New(&PyArray_Type, mdt->base.nfeat, dims, type_num, NULL,
                    data, 0, NPY_WRITEABLE, NULL);
  if (!obj) {
    g_free(dims);
    return NULL;
  }

  /* Ensure that the mdt.Table is kept around as long as the numpy object
     is alive. */
  Py_INCREF(mdt_pyobj);
  PyArray_BASE(obj) = mdt_pyobj;

  g_free(dims);
  return obj;
#else
  PyErr_SetString(PyExc_NotImplementedError,
                  "MDT was built without NumPy support");
  return NULL;
#endif
}
%}
