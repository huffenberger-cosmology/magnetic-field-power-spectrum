#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>

#include <numpy/ndarrayobject.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>

#include <lsstools.h>


static PyObject *Delta2k(PyObject *self, PyObject *args) {

  PyObject *N_array=NULL;
  PyObject *Delta_array=NULL;
  
  if (!PyArg_ParseTuple(args, "OO", &Delta_array, &N_array)) 
    return NULL;
  
  int *N = PyArray_DATA(N_array);
  double *Delta = PyArray_DATA(Delta_array);
  
  double *Deltak = calloc(3,sizeof(double));
  lsstools_Delta2k(N,Delta,Deltak);
  
  npy_intp dim = 3;
  PyObject *arr = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, Deltak);
  PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
   
  return(arr);
}



static PyMethodDef BpowspecMethods[] = {
  {"Delta2k",  Delta2k, METH_VARARGS,"Delta2k(double Delta[3], int N[3])"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef Bpowspec_module = {
    PyModuleDef_HEAD_INIT,
    "Bpowspec",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    BpowspecMethods
};


PyMODINIT_FUNC PyInit_Bpowspec(void)
{
  PyObject *m;

  m = PyModule_Create(&Bpowspec_module);
  
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!

  return(m);
}
