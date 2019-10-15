#include <stdio.h>
#include <stdlib.h>
#include <lsstools.h>

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

