#include <stdio.h>
#include <stdlib.h>
#include <lsstools.h>
#include <math.h>
#include <Bpowspec.h>


int lsstools_mapsize(int N[3]) {
  return(N[0]*N[1]*N[2]);
}


int lsstools_harmsize(int N[3]) {
  return(N[0]*N[1]*(N[2]/2+1));
}

void lsstools_size2Delta(double size[3], int N[3], double Delta[3]) {
  Delta[0] = size[0]/N[0];
  Delta[1] = size[1]/N[1];
  Delta[2] = size[2]/N[2];
}

void lsstools_Delta2k(int N[3],double Delta[3], double Deltak[3]) {
  // k = omega = 2 pi f
  // f = 1/N/Delta
  Deltak[0] = 2.0*M_PI/N[0]/Delta[0];
  Deltak[1] = 2.0*M_PI/N[1]/Delta[1];
  Deltak[2] = 2.0*M_PI/N[2]/Delta[2];
}  


void Bpowspec_build_projector(double P[9], double khat[3]){

  int i,j;

  for (j=0;j<3;j++)
    for (i=0;i<3;i++) {
      P[j*3+i] = (i==j) ? 1 : 0;
      P[j*3+i] -= khat[i]*khat[j];
      //printf("ij %d %d khat[i] khat[j] %f %f\n",i,j,khat[i],khat[j]);
    }

}


void Bpowspec_Pk2harm(gsl_spline *Pk, fftw_complex *harmx, fftw_complex *harmy, fftw_complex *harmz, int N[3], double kmax, double Deltak[3], gsl_rng *r) {
  // Based on a spline-interpolated power spectrum, return fourier space coefficients for the three components of the magnetic field
  // These are arrange in preparation for a c2r transform to a volume of size N[3].
  // N[] = Nz, Ny, Nx
  // Note that vectors ik, kvec, and khat (and thus the projection matrix Proj) also follow the zyx order.

  
  int ix,iy,iz,p,ik[3],idim,jdim;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  double kvec[3],khat[3],Proj[9];
  double k,normpower;
  int Nhalfx = N[2]/2+1;

  fftw_complex tmpharm[3],tmpBharm[3];
  
  double norm = pow(2*M_PI,3)/Deltak[0]/Deltak[1]/Deltak[2];

  for (iz=0;iz<N[0];iz++) {
    for (iy=0;iy<N[1];iy++) {
      for (ix=0;ix<Nhalfx;ix++) {
	p = iz*N[1]*Nhalfx+iy*Nhalfx+ix;
	
	ik[0] = (iz<N[0]/2) ? iz : iz-N[0];
	ik[1] = (iy<N[1]/2) ? iy : iy-N[1];
	ik[2] = ix;
	
	k = 0;
	for (idim = 0;idim < 3;idim++) {
	  kvec[idim] = Deltak[idim]*ik[idim];
	  k += kvec[idim]*kvec[idim];
	}
	k = sqrt(k);
	
	if (k>0.0) {
	  
	  for (idim = 0;idim < 3;idim++) {
	    khat[idim] = kvec[idim]/k;
	  }
	  
	  Bpowspec_build_projector(Proj,khat);
	} else {

	  for (idim = 0;idim < 9;idim++) {
	    Proj[idim] = 0.0;
	  }
	  for (idim = 0;idim < 3;idim++) {
	    Proj[idim*3+idim] = 1.0;
	  }
	  
	}
	
	//	printf("k = %e\n",k);

	
	normpower = ((k<=0.0)||(k > kmax)) ? 0.0 : gsl_spline_eval (Pk, k, acc)*norm;
	normpower = (normpower < 0.0) ? 0.0 : normpower;

	for (idim = 0;idim < 3;idim++) {
	  tmpharm[idim] =  sqrt(normpower)*(gsl_ran_gaussian(r,M_SQRT1_2) + I*gsl_ran_gaussian(r,M_SQRT1_2));
	}
	
	for (jdim = 0;jdim < 3;jdim++) {
	  tmpBharm[jdim] = 0.;
	  for (idim = 0;idim < 3;idim++) {
	    tmpBharm[jdim] += Proj[jdim*3+idim]*tmpharm[idim];
	  }
	}

	// Load into outputs, note the order zyx
	harmz[p] = tmpBharm[0];
	harmy[p] = tmpBharm[1];
	harmx[p] = tmpBharm[2];
	

	if ((iz==0)&&(iy==0)&&(ix==0)) {
	  printf("(%d %d %d)\n k=%f normpower=%f\n",iz,iy,ix,k,normpower);
	  printf("tmpharm[0] = (%f %f)\n",creal(tmpharm[0]),cimag(tmpharm[0]));
	  printf("tmpharm[1] = (%f %f)\n",creal(tmpharm[1]),cimag(tmpharm[1]));
	  printf("tmpharm[2] = (%f %f)\n",creal(tmpharm[2]),cimag(tmpharm[2]));
	  printf("tmpBharm[0] = (%f %f)\n",creal(tmpBharm[0]),cimag(tmpBharm[0]));
	  printf("tmpBharm[1] = (%f %f)\n",creal(tmpBharm[1]),cimag(tmpBharm[1]));
	  printf("tmpBharm[2] = (%f %f)\n",creal(tmpBharm[2]),cimag(tmpBharm[2]));
	}

	
      }
    }
  }
  harmx[0] = sqrt(2)*creal(harmx[0]);
  harmy[0] = sqrt(2)*creal(harmy[0]);
  harmz[0] = sqrt(2)*creal(harmz[0]);
  
  gsl_interp_accel_free(acc);
}


void lsstools_harm2map(fftw_complex *harm, double *map, int N[3], double Delta[3]) {
  int Nharm = lsstools_harmsize(N);
  fftw_complex *harmcopy = calloc(Nharm,sizeof(fftw_complex));
  memcpy(harmcopy,harm,Nharm*sizeof(fftw_complex));

  fftw_plan plan = fftw_plan_dft_c2r_3d(N[0], N[1], N[2], harmcopy, map, FFTW_ESTIMATE);
  
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  free(harmcopy);

  double Deltak[3];
  
  lsstools_Delta2k(N,Delta, Deltak);
  
  double Dvolk_2pi3 = Deltak[0]*Deltak[1]*Deltak[2]/pow(2*M_PI,3);
  double Npix = lsstools_mapsize(N);
  
  int p;
  
  for (p=0;p<Npix;p++) {
    map[p] *= Dvolk_2pi3;
  }
}
