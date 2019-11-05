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

static PyObject *harm2Pk(PyObject *self, PyObject *args) {

  PyObject *harmx_array=NULL;
  PyObject *harmy_array=NULL;
  PyObject *harmz_array=NULL;
  PyObject *Deltak_array=NULL;
  PyObject *k_array=NULL;


  if (!PyArg_ParseTuple(args, "OOOOO", &harmx_array, &harmy_array, &harmz_array, &Deltak_array, &k_array)) 
    return NULL;

  fftw_complex *harmx = PyArray_DATA(harmx_array);
  fftw_complex *harmy = PyArray_DATA(harmy_array);
  fftw_complex *harmz = PyArray_DATA(harmz_array);
  double *Deltak = PyArray_DATA(Deltak_array);
  double *k = PyArray_DATA(k_array);


  npy_intp *npyNharm = PyArray_DIMS(harmx_array);
  npy_intp *kbinsize = PyArray_DIMS(k_array);
  
  int nkbin = kbinsize[0] - 1;
  int Nharm[3] = {npyNharm[0],npyNharm[1],npyNharm[2]};
  int N[3] = {Nharm[0],Nharm[1],(Nharm[2]-1)*2};


  // printf("about to allocate.\n");

  gsl_histogram *Pkobs = gsl_histogram_alloc (nkbin);

  // printf("about to set ranges %d cl->n %d.\n", nlbin, cl->n);


  gsl_histogram_set_ranges (Pkobs, k, nkbin+1);

  // printf("about to evaluate cl.\n");
  Bpowspec_harm2Pk(harmx,harmy,harmz,Pkobs,N,Deltak);
 
  double *Pkout = calloc(nkbin,sizeof(double));


  bcopy(Pkobs->bin, Pkout, nkbin*sizeof(double));
  npy_intp npykbin[1] = {nkbin};


  PyObject *arr = PyArray_SimpleNewFromData(1,npykbin, NPY_DOUBLE, Pkout);
 
  PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
 
  gsl_histogram_free (Pkobs);

  return(arr);
}


static PyObject *scalarPk2harm(PyObject *self, PyObject *args) {

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
  fftw_complex *harm = calloc(harmsize,sizeof(fftw_complex));

  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set (r, seed);
  
  lsstools_Pk2harm(Pkspline, harm, N, kmax, Deltak, r);

  gsl_rng_free(r); 

  PyObject *arr = PyArray_SimpleNewFromData(3, Nharm, NPY_CDOUBLE, harm);

  PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);

  gsl_spline_free(Pkspline);

  return(arr);
}




static PyObject *kvecs(PyObject *self, PyObject *args) {

  PyObject *N_array=NULL;
  PyObject *Deltak_array=NULL;
  
  if (!PyArg_ParseTuple(args, "OO", &N_array, &Deltak_array)) 
    return NULL;

  double *Deltak = PyArray_DATA(Deltak_array);
  int *N = PyArray_DATA(N_array);

  printf("N = %d %d %d\n",N[0],N[1],N[2]);

  npy_intp Nharm[3] = {N[0],N[1],N[2]/2+1};
  int harmsize = lsstools_harmsize(N);
  double *kx = calloc(harmsize,sizeof(fftw_complex));
  double *ky = calloc(harmsize,sizeof(fftw_complex));
  double *kz = calloc(harmsize,sizeof(fftw_complex));

  Bpowspec_kvecs(kx, ky, kz, N, Deltak);

  PyObject *arrx = PyArray_SimpleNewFromData(3, Nharm, NPY_DOUBLE, kx);
  PyObject *arry = PyArray_SimpleNewFromData(3, Nharm, NPY_DOUBLE, ky);
  PyObject *arrz = PyArray_SimpleNewFromData(3, Nharm, NPY_DOUBLE, kz);

  PyArray_ENABLEFLAGS((PyArrayObject *)arrx, NPY_OWNDATA);
  PyArray_ENABLEFLAGS((PyArrayObject *)arry, NPY_OWNDATA);
  PyArray_ENABLEFLAGS((PyArrayObject *)arrz, NPY_OWNDATA);

  PyObject *kvecstuple = PyTuple_Pack(3,arrx,arry,arrz);
  
  return(kvecstuple);
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

static PyObject *map2harm(PyObject *self, PyObject *args) {

  PyObject *map_array=NULL;
  PyObject *Delta_array=NULL;

  if (!PyArg_ParseTuple(args, "OO", &map_array, &Delta_array)) 
    return NULL;

  double *map = PyArray_DATA(map_array);

  npy_intp *npyN = PyArray_DIMS(map_array);
  int N[3] = {npyN[0],npyN[1],npyN[2]};

  double *Delta = PyArray_DATA(Delta_array);

  npy_intp Nharm[3] = {N[0],N[1],N[2]/2+1};
  fftw_complex *harm = calloc(Nharm[0]*Nharm[1]*Nharm[2],sizeof(fftw_complex));
  
  lsstools_map2harm(map,harm,N,Delta);

  PyObject *arr = PyArray_SimpleNewFromData(3, Nharm, NPY_CDOUBLE, harm);

  PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);

  
  return(arr);
}


static PyMethodDef BpowspecMethods[] = {
  {"Delta2k",  Delta2k, METH_VARARGS,"Delta2k(double Delta[3], int N[3])"},
  {"build_projector", build_projector, METH_VARARGS,"build_projector(double khat[3])"},
  {"Pk2harm",  Pk2harm, METH_VARARGS,"Pk2harm(k[], Pk[], N[3] (int32), kmax, Deltak, random seed (int))"},
  {"scalarPk2harm",  scalarPk2harm, METH_VARARGS,"scalarPk2harm(k[], Pk[], N[3] (int32), kmax, Deltak, random seed (int))"},
  {"harm2Pk",  harm2Pk, METH_VARARGS,"harm2Pk(harmx,harmy,harmz (complex), Deltak[2], kbin_edges[] (double))"},
  //  {"scalarharm2Pk",  scalarharm2Pk, METH_VARARGS,"scalarharm2Pk(harmx (complex), Deltak[2], kbin_edges[] (double))"},
  {"harm2map", harm2map, METH_VARARGS,"harm2map(harm[] (complex), Delta[3])"},
  {"map2harm", map2harm, METH_VARARGS,"map2harm(map[] (double), Delta[3])"},
  {"kvecs", kvecs, METH_VARARGS,"kvecs(N[3] (int32), Deltak)"},
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
