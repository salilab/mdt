%{
/* Get a user-defined property */
static float *python_cb_get_property(gpointer data,
                                     const struct mod_alignment *aln, int iseq,
                                     const struct mdt_library *mlib,
                                     const struct mod_libraries *libs,
                                     GError **err)
{
  PyObject *func, *arglist, *result;

  func = data;
  arglist = Py_BuildValue("(OiOO)", aln->scriptobj, iseq, mlib->scriptobj,
                          libs->scriptobj);
  if (!arglist) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "Python error");
    return NULL;
  }

  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);

  if (result) {
    int reslen;
    float *res = to_list_float(result, 0, &reslen, "result");
    Py_DECREF(result);
    if (!res) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "Python error");
    }
    return res;
  } else {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "Python error");
    return NULL;
  }
}

/* Checks that the provided Python object is a function, and sets
   parameters necessary to call it from a C callback function. The
   Python reference count is increased to prevent the function from
   being garbage collected. */
static gboolean set_python_callback(void **funcpt, gpointer *data,
                                    PyObject *pyfunc, void *cfunc)
{
  *funcpt = cfunc;
  *data = pyfunc;
  if (!PyCallable_Check(pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return FALSE;
  }
  Py_INCREF(pyfunc);
  return TRUE;
}

/* Helper function to release the reference to a Python function; used
   by set_python_callback_free() */
static void python_callback_decref(gpointer data)
{
  PyObject *obj = data;
  Py_DECREF(obj);
}

/* Set up a Python callback, similar to set_python_callback(). Additionally,
   a 'free' function is registered to release the reference to the Python
   function when it is no longer required. */
static gboolean set_python_callback_free(void **funcpt, gpointer *data,
                                         GDestroyNotify *freept,
                                         PyObject *pyfunc, void *cfunc)
{
  if (set_python_callback(funcpt, data, pyfunc, cfunc)) {
    *freept = python_callback_decref;
    return TRUE;
  } else {
    return FALSE;
  }
}
%}

/* User-defined property */
%typemap(in) (mdt_cb_get_property get_property, gpointer data,
              GDestroyNotify freefunc) {
  if (!set_python_callback_free((void **)&$1, &$2, &$3, $input,
                                python_cb_get_property)) {
    SWIG_fail;
  }
}
