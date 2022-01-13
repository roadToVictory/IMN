#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdio>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>
#include "mgmres.h"
#include "mgmres.c"


#define delta 0.1
#define itr_max 500
#define mr 500
#define tol_abs 1e-8
#define tol_rel 1e-8

double ro1(int nx, int ny, int i, int j, double sgm){
  return exp(  -(pow((i*delta-0.25*nx*delta),2))/(pow(sgm,2))  -(pow((j*delta-0.5*ny*delta),2))/(pow(sgm,2)));
}

double ro2(int nx, int ny, int i, int j, double sgm){
  return -exp(  -(pow((i*delta-0.75*nx*delta),2))/(pow(sgm,2))  -(pow((j*delta-0.5*ny*delta),2))/(pow(sgm,2)));
}

void Poisson(int nx, int ny, double eps1, double eps2, double V1, double V2, double V3, double V4, FILE *fa, FILE* fb, bool zad6){
  FILE *V_=fopen(("V"+std::to_string(nx)+"-"+std::to_string(zad6)+"-"+std::to_string(eps2)+".dat").c_str(), "w");
  if(!V_) return;

  int N = (nx+1)*(ny+1);
  double* a = new double[5*N];
  int* ja = new int[5*N];
  int* ia = new int[N+1];
  double* b = new double[N];
  double* V = new double[N];

  for(int i=0; i<=N; i++)
    ia[i]= -1;
  
  int k = -1; 
  int i, j, brzeg=0, nz_num;
  double sgm = delta*nx/10.;
  double* E = new double[N+nx+1];

  for(int l=0; l<N+nx+1; l++){
    j=floor(l/(nx+1));          //
    i=l-j*(nx+1);             //

    E[l] = (i<=nx/2)?eps1:eps2;
  }

  for(int l=0; l<N; l++){
    j=floor(l/(nx+1));          //
    i=l-j*(nx+1);    
    brzeg=0;
    double vb =0.; //potencjał na brzegu

    if(i==0){ //lewy brzeg
      brzeg = 1; 
      vb=V1;
    }

    if(j==ny){ //gorny brzeg
      brzeg = 1; 
      vb=V2;
    }

    if(i==nx){ //prawy brzeg
      brzeg = 1; 
      vb=V3;
    }

    if(j==0){ //dolny brzeg
      brzeg = 1; 
      vb=V4;
    }

    if(zad6)
      b[l] = -(ro1(nx, ny, i, j, sgm) + ro2(nx, ny, i, j, sgm));
    else b[l]=0.;

    if(brzeg==1)
      b[l] = vb;

    ia[l] = -1;

    if(l-nx-1>=0 && brzeg==0){  //lewa skrajna
      k++;

      if(ia[l]<0)
        ia[l]=k;

      a[k] = E[l]/pow(delta,2);
      ja[k] = l-nx-1;
    }

    if(l-1>=0 && brzeg==0){  //poddiagonala
      k++;

      if(ia[l]<0)
        ia[l]=k;

      a[k] = E[l]/pow(delta,2);
      ja[k] = l-1;
    }

    //diagonala
    k++;
    if(ia[l]<0)
      ia[l]=k;
    if(brzeg==0)
      a[k] = -(2*E[l]+E[l+1]+E[l+nx+1])/(pow(delta,2));
    else a[k]=1;

    ja[k]=l;

    if(l<N && brzeg==0){ //naddiagonala
      k++;
      a[k] = E[l+1]/pow(delta,2);
      ja[k] = l+1;
    }

    if(l<(N-nx-1) && brzeg==0){ //prawa skrajna
      k++;
      a[k] = E[l+nx+1]/pow(delta,2);
      ja[k]= l+nx+1;
    }
  }

  nz_num=k+1;
  ia[N] = nz_num;

  pmgmres_ilu_cr(N, nz_num, ia, ja, a,V, b, itr_max, mr, tol_abs, tol_rel);

  for(int l=0; l<N; l++ ){
    j = floor( l / (nx + 1) );
    i = l - j * (nx + 1);
    fprintf(V_, "%d\t%g\t%g\t%g\n", l, i*delta, j*delta, V[l]);
    fprintf(fa, "%d\t%d\t%d\t%g\n", l, i, j, a[l]);
    fprintf(fb, "%d\t%d\t%d\t%g\n", l, i, j, b[l]);

    if(i==nx)
      fprintf(V_, "\n");
  }

  fclose(V_);
  delete[] a; delete[] ja; delete[] ia; delete[] b; delete[] V;
}


int main(){
//warunki początkowe
  int nx =4, ny=4;
  double eps1 = 1., eps2=1.;
  double V1 = 10., V3 = 10., V2=-10., V4=-10.;

  FILE *fa = fopen("A.dat", "w");
  FILE *fb = fopen("b.dat", "w");
  if(!fa || !fb) return 0;

  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, false);

  nx = ny =50;
  std::cout<<nx<<" "<<ny;
  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, false);

  nx = ny =100;
  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, false);

  nx = ny =200;
  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, false);

  //rozklad
  nx=ny=100;
  V1=V2=V3=V4=0.;

  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, true);

  eps1=1; eps2=2;
  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, true);

  eps1=1; eps2=10;
  Poisson(nx, ny, eps1, eps2, V1, V2, V3, V4, fa, fb, true);

  fclose(fa); fclose(fb);
  return 0;
}