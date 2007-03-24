%{
#ifdef SWIGPYTHON
/* Generic error */
static PyObject *mdterror;

/** Raise an exception if an error code was returned. */
static void handle_error(GError *err)
{
  PyObject *py_err_type = mdterror;
  if (err->domain == MDT_ERROR) {
    switch(err->code) {
    case MDT_ERROR_IO:
      py_err_type = PyExc_IOError;
      break;
    case MDT_ERROR_VALUE:
      py_err_type = PyExc_ValueError;
      break;
    case MDT_ERROR_INDEX:
      py_err_type = PyExc_IndexError;
      break;
    }
  }
  PyErr_SetString(py_err_type, err->message);
  g_error_free(err);
}
#endif /* SWIGPYTHON */
%}

%typemap(in, numinputs=0) GError **err (GError *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout) GError **err {
  if (*$1) {
    handle_error(*$1);
#ifdef SWIGPYTHON
    Py_DECREF(resultobj);
#endif
    SWIG_fail;
  }
}
