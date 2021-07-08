/*
    This code is part of magnetic-field-power-spectrum.
    Copyright (C) 2021 Kevin M. Huffenberger, khuffenberger@fsu.edu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <fftw3.h>
#include <gsl/gsl_histogram.h>

void Bpowspec_build_projector(double P[9], double khat[3]);
void Bpowspec_Pk2harm(gsl_spline *Pk, fftw_complex *harmx, fftw_complex *harmy, fftw_complex *harmz, int N[3], double kmax, double Deltak[3], gsl_rng *r);
void Bpowspec_harm2Pk(fftw_complex *harmx, fftw_complex *harmy, fftw_complex *harmz, gsl_histogram *Pkobs, int N[3], double Deltak[3]);



void Bpowspec_kvecs(double *kx, double *ky, double *kz, int N[3], double Deltak[3]);
