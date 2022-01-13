#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdio>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>

//parametry
#define delta 0.01
#define p 1.
#define um 1.
#define nx 200
#define ny 90
#define i1 50
#define j1 55
#define IT_MAX 20000

void WB_psi(double Qwe, double Qwy, double* y, double** Psi){
  for(int j=j1; j<=ny;j++)     //brzeg A
    Psi[0][j]= Qwe/(2*um) * (pow(y[j],3)/3 - pow(y[j],2)/2 *(y[j1]+y[ny])+ y[j]*y[j1]*y[ny]);

  for(int j=0; j<=ny; j++)      //brzeg C
    Psi[nx][j] = Qwy/(2*um) * (pow(y[j],3)/3 - (pow(y[j],2)/2)* y[ny] ) + (Qwe*pow(y[j1],2)*(-y[j1]+3*y[ny]))/(12*um);
  
  for(int i=1; i<=nx-1; i++)    //brzeg B
    Psi[i][ny] = Psi[0][ny];

  for(int i=i1; i<=nx-1; i++)    //brzeg D
    Psi[i][0] = Psi[0][j1];

  for(int j=1; j<=j1;j++)    //brzeg E
    Psi[i1][j] = Psi[0][j1];

  for(int i=1; i<=i1; i++)    //brzeg F
    Psi[i][j1] = Psi[0][j1];
}

void WB_GZeta(double Qwe, double Qwy, double* y, double** GZeta, double **psi){
  for(int j=j1; j<=ny; j++) //brzeg A
    GZeta[0][j] = Qwe/(2*um)*(2*y[j]-y[j1]-y[ny]);

  for(int j=0; j<=ny; j++) //brzeg C
    GZeta[nx][j] = Qwy/(2*um)*(2*y[j]-y[ny]);
  
  for(int i=1; i<=nx-1; i++)    //brzeg B
    GZeta[i][ny] = 2/pow(delta,2) * (GZeta[i][ny-1]- GZeta[i][ny]);

  for(int i=i1+1; i<=nx-1; i++)    //brzeg D
    GZeta[i][0] = 2/pow(delta,2) * (GZeta[i][1]- GZeta[i][0]);

  for(int j=1; j<=j1-1;j++)    //brzeg E
    GZeta[i1][j] = 2/pow(delta,2) * (GZeta[i1+1][j]- GZeta[i1][j]);

  for(int i=1; i<=i1; i++)    //brzeg F
    GZeta[i][j1] =( 2/pow(delta,2)) * (GZeta[i][j1+1]- GZeta[i][j1]);

  GZeta[i1][j1] = (GZeta[i1-1][j1]+GZeta[i1][j1-1])/2;
}

bool czyNaBrzegu(int i, int j){
    return (i!=0 && j !=0 && i!=nx && j!=ny && !(i<=i1 && j==j1) && !(i==i1 && j<=j1) && !(i<=i1 && j < j1));
}

void Relaksacja(double Qwe, FILE* file){

  double** Psi = new double*[nx+1]; 
  double** GZeta = new double*[nx+1]; 
  double** u = new double*[nx+1]; 
  double** v = new double*[nx+1]; 
  double omega;

  for(int i=0; i<=nx;i++){
    Psi[i] = new double[ny+1]; 
    GZeta[i] = new double[ny+1]; 
    u[i] = new double[ny+1]; 
    v[i] = new double[ny+1]; 
  }

  for(int i=0; i<=nx;i++)
    for(int j=0; j<=ny;j++){
      Psi[i][j]=0.;
      GZeta[i][j]=0.;
      u[i][j]=0.;
      v[i][j]=0.;
    }

  double Qwy = Qwe*(pow((delta*ny),3) - pow((delta*j1),3) -3*delta*j1*pow((delta*ny),2) + 3*pow((delta*j1),2)*delta*ny )/(pow((delta*ny),3));

  double* y = new double [ny + 1];
  for (int i = 0; i <= ny; i++)
        y[i] =  delta * i;
	

  WB_psi(Qwe, Qwy, y, Psi);

  for(int it=1; it<=IT_MAX; it++){//
    omega = it<2000? 0:1;

    for(int i=1; i<=nx-1; i++)
      for(int j=1; j<=ny-1; j++)
          if(czyNaBrzegu(i, j)){
            Psi[i][j] = (Psi[i+1][j]+Psi[i-1][j]+Psi[i][j+1]+Psi[i][j-1]- pow(delta,2)*GZeta[i][j])/4.;
            GZeta[i][j] = (GZeta[i+1][j]+GZeta[i-1][j]+GZeta[i][j+1]+GZeta[i][j-1])/4. - omega*p/(16*um)*( (Psi[i][j+1]-Psi[i][j-1])*(GZeta[i+1][j]-GZeta[i-1][j]) - (Psi[i+1][j]-Psi[i-1][j]) * (GZeta[i][j+1]-GZeta[i][j-1]) );
//pochodne zastępujemy ilorazami roznicowymi
            u[i][j] = (Psi[i][j+1] - Psi[i][j-1]) / (2*delta);
            v[i][j] = -(Psi[i+1][j] - Psi[i-1][j]) / (2*delta);
          }
    WB_GZeta(Qwe, Qwy, y, GZeta, Psi);

    double er=0.;
    for(int i=1; i<=nx-1;i++)
      er+= Psi[i+1][j1+2]+ Psi[i-1][j1+2] + Psi[i][j1+3] + Psi[i][j1+1] - 4*Psi[i][j1+2]-pow(delta,2)*GZeta[i][j1+2];
  }

  for (int i = 1; i<=nx-1; i++){
    for (int j = 1; j<=ny-1; j++)
      if (i <= i1 && j <= j1){
               fprintf(file, "%f\t %f\t %g\t %g\t %g\t %g\n", i * delta, j * delta, Psi[i1][j1], GZeta[i1][j1], u[i][j], v[i][j]);
            }
			else
                fprintf(file, "%f\t %f\t %g\t %g\t %g\t %g\n", i * delta, j * delta, Psi[i][j], GZeta[i][j], u[i][j], v[i][j]);

        fprintf(file, "\n");
  }

  for(int i=0; i<=nx;i++){
      delete[] Psi[i]; 
      delete[] GZeta[i]; 
      delete[] u[i];
      delete[] v[i];
    }
  delete[] GZeta; delete[] Psi; delete[] u; delete[] v; 
}


int main(){
  double Q = -1000.;
  FILE *f1= fopen(("Q_"+std::to_string(Q)+".dat").c_str(), "w");
  if(f1) Relaksacja(Q, f1);

  Q = -4000.;
  FILE *f2= fopen(("Q_"+std::to_string(Q)+".dat").c_str(), "w");
  if(f2) Relaksacja(Q, f2);

  Q = 4000.;
  FILE *f3= fopen(("Q_"+std::to_string(Q)+".dat").c_str(), "w");
  if(f3) Relaksacja(Q, f3);
  
  fclose(f1); fclose(f2); fclose(f3);
  return 0;
}