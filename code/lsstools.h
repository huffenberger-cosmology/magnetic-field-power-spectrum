#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <complex.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <fftw3.h>
#include <string.h>

int lsstools_mapsize(int N[3]);
int lsstools_harmsize(int N[3]);

void lsstools_size2Delta(double size[3], int N[3], double Delta[3]);
void lsstools_Delta2k(int N[3],double Delta[3], double Deltak[3]);



void lsstools_Pk2harm(gsl_spline *Pk, fftw_complex *harm, int N[3], double kmax, double Deltak[3], gsl_rng *r);
void lsstools_harm2Pk(fftw_complex *harm, gsl_histogram *Pk, int N[3], double Deltak[3]);



void lsstools_harm2map(fftw_complex *harm, double *map, int N[3], double Delta[3]);
void lsstools_map2harm( double *map, fftw_complex *harm, int N[3], double Delta[3]);

void lsstools_harm_filter(gsl_spline *fk, fftw_complex *harm, int N[3], double kmax, double Deltak[3]);

/*
void lsstools_discrete_kernel(gsl_spline *kernel, double *F, int N, double Delta, double offset);
void lsstools_integratex(gsl_spline *kernel, double *map, int N[3], double Delta[3],double offsetx, double *map2d);
void lsstools_projectx(gsl_spline *kernel, double *map, int N[3], double Delta[3],double offsetx, double *map2d);


 int lsstools_map2fits(char *filename, double *map, int N[3]);
int lsstools_fits2map(char *filename, double **map, int N[3]);
*/
