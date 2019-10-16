#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>

#include <numpy/ndarrayobject.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>



#include <Bpowspec.h>
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

static PyObject *build_projector(PyObject *self, PyObject *args) {

  PyObject *khat_array=NULL;
  
  if (!PyArg_ParseTuple(args, "O", &khat_array)) 
    return NULL;
 
  double *P = calloc(9,sizeof(double));

  double *khat = PyArray_DATA(khat_array);

  
  Bpowspec_build_projector(P, khat);

  int dim = 2;
  npy_intp sz[2] = {3,3};
  PyObject *arr = PyArray_SimpleNewFromData(dim, sz, NPY_DOUBLE, P);
  PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);

  return(arr);
}


static PyObject *Pk2harm(PyObject *self, PyObject *args) {

  PyObject *k_array=NULL;
  PyObject *Pk_array=NULL;
  PyObject *N_array=NULL;
  PyObject *Deltak_array=NULL;
  double kmax;
  int seed;
  
  
  if (!PyArg_ParseTuple(args, "OOOdOi", &k_array, &Pk_array, &N_array, &kmax, &Deltak_array, &seed)) 
    return NULL;

  double *k = PyArray_DATA(k_array);
  double *Pk = PyArray_DATA(Pk_array);
  
  npy_intp *Nkptr = PyArray_DIMS(k_array);
  int Nk = (int) *Nkptr;
  double *Deltak = PyArray_DATA(Deltak_array);

  gsl_spline *Pkspline = gsl_spline_alloc (gsl_interp_linear, Nk);

  
  gsl_spline_init (Pkspline, k, Pk, Nk);

  int *N = PyArray_DATA(N_array);

  printf("N = %d %d %d\n",N[0],N[1],N[2]);

  npy_intp Nharm[3] = {N[0],N[1],N[2]/2+1};
  int harmsize = lsstools_harmsize(N);
  fftw_complex *harmx = calloc(harmsize,sizeof(fftw_complex));
  fftw_complex *harmy = calloc(harmsize,sizeof(fftw_complex));
  fftw_complex *harmz = calloc(harmsize,sizeof(fftw_complex));

  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set (r, seed);
  
  Bpowspec_Pk2harm(Pkspline, harmx, harmy, harmz, N, kmax, Deltak, r);

  gsl_rng_free(r); 

  PyObject *arrx = PyArray_SimpleNewFromData(3, Nharm, NPY_CDOUBLE, harmx);
  PyObject *arry = PyArray_SimpleNewFromData(3, Nharm, NPY_CDOUBLE, harmy);
  PyObject *arrz = PyArray_SimpleNewFromData(3, Nharm, NPY_CDOUBLE, harmz);

  PyArray_ENABLEFLAGS((PyArrayObject *)arrx, NPY_OWNDATA);
  PyArray_ENABLEFLAGS((PyArrayObject *)arry, NPY_OWNDATA);
  PyArray_ENABLEFLAGS((PyArrayObject *)arrz, NPY_OWNDATA);

  gsl_spline_free(Pkspline);

  PyObject *harmtuple = PyTuple_Pack(3,arrx,arry,arrz);
  
  return(harmtuple);
}

static PyObject *harm2map(PyObject *self, PyObject *args) {

  PyObject *harm_array=NULL;
  PyObject *Delta_array=NULL;

  if (!PyArg_ParseTuple(args, "OO", &harm_array, &Delta_array)) 
    return NULL;
  
  fftw_complex *harm = PyArray_DATA(harm_array);

  npy_intp *Nharm = PyArray_DIMS(harm_array);

  printf("harm2map harmdim %d %d %d\n",(int)Nharm[0],(int)Nharm[1],(int)Nharm[2]);

  int N[3] = {Nharm[0],Nharm[1],2*(Nharm[2]-1)};

  double *Delta = PyArray_DATA(Delta_array);

  double *map = calloc(N[0]*N[1]*N[2],sizeof(double));
  
  lsstools_harm2map(harm, map, N, Delta);
    
  npy_intp npyN[] = {N[0],N[1],N[2]};

  PyObject *arr = PyArray_SimpleNewFromData(3, npyN, NPY_DOUBLE, map);
  PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);

  return(arr);

}



static PyMethodDef BpowspecMethods[] = {
  {"Delta2k",  Delta2k, METH_VARARGS,"Delta2k(double Delta[3], int N[3])"},
  {"build_projector", build_projector, METH_VARARGS,"build_projector(double khat[3])"},
  {"Pk2harm",  Pk2harm, METH_VARARGS,"Pk2harm(k[], Pk[], N[3] (int32), kmax, Deltak, random seed (int))"},
  {"harm2map", harm2map, METH_VARARGS,"harm2map(harm[] (complex), Delta[3])"},
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
