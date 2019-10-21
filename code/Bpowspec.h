#include <fftw3.h>

void Bpowspec_build_projector(double P[9], double khat[3]);
void Bpowspec_Pk2harm(gsl_spline *Pk, fftw_complex *harmx, fftw_complex *harmy, fftw_complex *harmz, int N[3], double kmax, double Deltak[3], gsl_rng *r);
void Bpowspec_kvecs(double *kx, double *ky, double *kz, int N[3], double Deltak[3]);
