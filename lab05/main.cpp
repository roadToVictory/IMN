#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdio>
#include <cmath>
#include <math.h>
#include <string>

//parametry zadania
double delta = 0.2;
int nx=128;
int ny=128;
double xmax = delta*nx;
double ymax = delta*ny;
int iter=0;

void WB(double **V){ //warunki brzegowe
  for(int i=0; i<nx+1; i++){
    V[0][i] = V[nx][i] = std::sin(M_PI*delta*i/ymax ); //VB_1, VB_3
    V[i][ny] = -std::sin(2*M_PI*delta*i/xmax ); //VB_2
    V[i][0] = std::sin(2*M_PI*i*delta/xmax); //VB_4
  }
}

void Poisson(){
  double** V = new double* [nx+1];
  double S_it = 0.;
  double S_it_p = 0.;
  int k[]={16,8,4,2,1};

  for(int i=0; i<nx+1; i++){
    V[i] = new double[ny+1];
    for(int j=0; j<ny+1; j++)
      V[i][j]=0.;
  }
  WB(V);

  for(int a=0; a<5; a++){
    int k1 = k[a];
    FILE* Va = fopen(("V"+std::to_string(k1)+".dat").c_str(), "w");
    FILE* Sa = fopen(("S"+std::to_string(k1)+".dat").c_str(), "w");

    for(int i=0; i<=nx-k1; i+=k1)
      for(int j=0; j<=nx-k1; j+=k1)
        S_it += (pow(k1*delta,2)/2) * (  pow(  ((V[i+k1][j]-V[i][j])/(2*k1*delta)) + ((V[i+k1][j+k1]-V[i][j+k1])/(2*k1*delta)) ,2) + pow(  ((V[i][j+k1]-V[i][j])/(2*k1*delta)) + ((V[i+k1][j+k1]-V[i+k1][j])/(2*k1*delta)) ,2));

     while(1){
        iter++;

        for(int i=k1; i<=nx-k1; i+=k1)
          for(int j=k1; j<=ny-k1; j+=k1)
            V[i][j] = 0.25*(V[i+k1][j]+V[i-k1][j]+V[i][j+k1]+V[i][j-k1]);

        S_it_p=S_it;
        S_it=0.;
        for(int i=0; i<=nx-k1; i+=k1)
          for(int j=0; j<=nx-k1; j+=k1)
            S_it += (pow(k1*delta,2)/2) * (  pow(  ((V[i+k1][j]-V[i][j])/(2*k1*delta)) + ((V[i+k1][j+k1]-V[i][j+k1])/(2*k1*delta)) ,2) + pow(  ((V[i][j+k1]-V[i][j])/(2*k1*delta)) + ((V[i+k1][j+k1]-V[i+k1][j])/(2*k1*delta)) ,2));

        fprintf(Sa, "%d\t%g\n", iter, S_it);
		
        if(fabs((S_it-S_it_p)/S_it_p)<1e-8) break;
      }

      for(int i=0; i<nx+1; i+=k1){
        for(int j=0; j<ny+1; j+=k1){
          fprintf(Va,"%d\t%d\t%g\n", i,j,V[i][j]);
        }
        fprintf(Va, "\n");
      }

  if(k1!=1){
    for(int i=0; i<=nx-k1; i+=k1){
        for(int j=0; j<=ny-k1; j+=k1){
          V[i+k1/2][j+k1/2] = 0.25*(V[i][j]+V[i+k1][j]+V[i][j+k1]+V[i+k1][j+k1]); //środkowy

          if(i+k1 != nx) 
            V[i+k1][j+k1/2] = 0.5*(V[i+k1][j]+V[i+k1][j+k1]); //prawy

          if(j+k1 != ny)
            V[i+k1/2][j+k1] = 0.5*(V[i][j+k1]+V[i+k1][j+k1]); //gorny

          if(j != 0)
            V[i+k1/2][j] = 0.5*(V[i][j]+V[i+k1][j]); //dolny

          if(i != 0)
            V[i][j+k1/2] = 0.5*(V[i][j]+V[i][j+k1]); //lewy
        }
    }
  }

    fclose(Va); fclose(Sa); 
  }

  for(int i=0; i<nx+1; i++) delete[] V[i];
  delete[] V;
}

int main(){
  Poisson();
  return 0;
}