#ifdef SWIGPYTHON
%define TO_LIST(type, checkfn, convertfn, errmsg)
%{
static type *to_list_ ## type ##(PyObject *pyinput, int fixsize, int *sizevar, char *displayname) {
  Py_ssize_t seqlen;
  int i, intseqlen;
  type *outlist;
  /* Treat a single value the same as a 1-element list */
  if (checkfn(pyinput) && sizevar) {
    outlist = g_malloc(sizeof(type));
    *sizevar = 1;
    outlist[0] = (type)convertfn(pyinput);
    return outlist;
  }
#if PY_VERSION_HEX < 0x03000000
  if (!PySequence_Check(pyinput) || PyString_Check(pyinput)) {
#else
  if (!PySequence_Check(pyinput) || PyUnicode_Check(pyinput)
      || PyBytes_Check(pyinput)) {
#endif
    PyErr_Format(PyExc_ValueError, "%s should be a sequence", displayname);
    return NULL;
  }
  seqlen = PySequence_Length(pyinput);
  /* Make sure length fits into an int (necessary for correct functioning
     with Python 2.5 on 64-bit platforms) */
  if (seqlen > INT_MAX) {
    PyErr_Format(PyExc_ValueError, "sequence %s length exceeds maximum",
                 displayname);
    return NULL;
  }
  intseqlen = (int)seqlen;
  if (sizevar) {
    *sizevar = intseqlen;
  } else if (intseqlen != fixsize) {
    PyErr_Format(PyExc_ValueError, "%s must be a sequence of length %d",
                 displayname, fixsize);
    return NULL;
  }
  /* malloc(0) is undefined, so make sure we use a non-zero size */
  outlist = g_malloc(sizeof(type) * (intseqlen == 0 ? 1 : intseqlen));

  for (i = 0; i < intseqlen; i++) {
    PyObject *o = PySequence_GetItem(pyinput, i);
    if (checkfn(o)) {
      outlist[i] = (type)convertfn(o);
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      PyErr_Format(PyExc_ValueError, errmsg, displayname, i);
      g_free(outlist);
      return NULL;
    }
  }
  return outlist;
}
%}
%enddef

TO_LIST(float, PyNumber_Check, PyFloat_AsDouble, "%s[%d] should be a number")
TO_LIST(int, PyInt_Check, PyInt_AsLong, "%s[%d] should be an integer")


%typemap(in) (const float VARLIST[], int N_VARLIST) {
    $1 = to_list_float($input, 0, &$2, "$1_name");
    if (!$1) SWIG_fail;
}

%typemap(in) (const int VARLIST[], int N_VARLIST) {
    $1 = to_list_int($input, 0, &$2, "$1_name");
    if (!$1) SWIG_fail;
}

%typemap(in) const int [ANY] {
    $1 = to_list_int($input, $1_dim0, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
#endif 

%typemap(freearg) (const float VARLIST[], int N_VARLIST) {
    g_free($1);
}
%typemap(arginit) (const float VARLIST[], int N_VARLIST) {
  $1 = NULL;
}
%typemap(freearg) (const int VARLIST[], int N_VARLIST) {
    g_free($1);
}
%typemap(arginit) (const int VARLIST[], int N_VARLIST) {
  $1 = NULL;
}
%typemap(freearg) const int [ANY] {
    g_free($1);
}
%typemap(arginit) const int [ANY] {
  $1 = NULL;
}
