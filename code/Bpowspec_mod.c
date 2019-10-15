#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>

static PyMethodDef BpowspecMethods[] = {
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initBpowspec(void)
{
  (void) Py_InitModule("Bpowspec", BpowspecMethods);
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
}
