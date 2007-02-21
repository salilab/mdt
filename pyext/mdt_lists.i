#ifdef SWIGPYTHON
%define TO_VARLIST(type, checkfn, convertfn, errmsg)
%inline %{
static type *to_varlist_ ## type ##(PyObject *pyinput, int *sizevar, char *displayname) {
  Py_ssize_t seqlen;
  int i, intseqlen;
  type *outlist;
  /* Treat a single value the same as a 1-element list */
  if (checkfn(pyinput)) {
    outlist = malloc(sizeof(type));
    *sizevar = 1;
    outlist[0] = (type)convertfn(pyinput);
    return outlist;
  }
  if (!PySequence_Check(pyinput) || PyString_Check(pyinput)) {
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
  *sizevar = intseqlen = (int)seqlen;
  /* malloc(0) is undefined, so make sure we use a non-zero size */
  outlist = malloc(sizeof(type) * (intseqlen == 0 ? 1 : intseqlen));

  for (i = 0; i < intseqlen; i++) {
    PyObject *o = PySequence_GetItem(pyinput, i);
    if (checkfn(o)) {
      outlist[i] = (type)convertfn(o);
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      PyErr_Format(PyExc_ValueError, errmsg, displayname, i);
      free(outlist);
      return NULL;
    }
  }
  return outlist;
}
%}
%enddef

TO_VARLIST(float, PyNumber_Check, PyFloat_AsDouble, "%s[%d] should be a number")
TO_VARLIST(int, PyInt_Check, PyInt_AsLong, "%s[%d] should be an integer")


%typemap(in) (const float VARLIST[], int N_VARLIST) {
    $1 = to_varlist_float($input, &$2, "$1_name");
    if (!$1) SWIG_fail;
}

%typemap(in) (const int VARLIST[], int N_VARLIST) {
    $1 = to_varlist_int($input, &$2, "$1_name");
    if (!$1) SWIG_fail;
}
#endif 

%typemap(freearg) (const float VARLIST[], int N_VARLIST) {
    if ($1) free($1);
}
%typemap(arginit) (const float VARLIST[], int N_VARLIST) {
  $1 = NULL;
}
%typemap(freearg) (const int VARLIST[], int N_VARLIST) {
    if ($1) free($1);
}
%typemap(arginit) (const int VARLIST[], int N_VARLIST) {
  $1 = NULL;
}
