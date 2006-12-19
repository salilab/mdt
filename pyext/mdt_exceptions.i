%{
#ifdef SWIGPYTHON
/* Generic error */
static PyObject *mdterror;


/** Raise an exception if an error code was returned. */
static int handle_error(void)
{
  char *mod_err = get_mod_error();

  if (mod_err) {
    PyObject *py_err_type;
    int mod_err_type = get_mod_error_type();
    switch(mod_err_type) {
    case ME_EOF:
      py_err_type = PyExc_EOFError;
      break;
    case ME_INDEX:
      py_err_type = PyExc_IndexError;
      break;
    case ME_IO:
      py_err_type = PyExc_IOError;
      break;
    case ME_MEMORY:
      py_err_type = PyExc_MemoryError;
      break;
    case ME_NOTIMP:
      py_err_type = PyExc_NotImplementedError;
      break;
    case ME_TYPE:
      py_err_type = PyExc_TypeError;
      break;
    case ME_VALUE:
      py_err_type = PyExc_ValueError;
      break;
    case ME_ZERODIV:
      py_err_type = PyExc_ZeroDivisionError;
      break;
    default:
      py_err_type = mdterror;
    }
    PyErr_SetString(py_err_type, mod_err);
  } else if (!PyErr_Occurred()) {
    PyErr_SetString(mdterror,
                    "INTERNAL ERROR: error code set, but no error information");
  }
}
#endif /* SWIGPYTHON */
%}

%typemap(in, numinputs=0) int *ierr (int temp) {
  $1 = &temp;
}

%typemap(argout) int *ierr {
  if (*$1 != 0) {
    handle_error();
#ifdef SWIGPYTHON
    Py_DECREF(resultobj);
#endif
    SWIG_fail;
  }
}
