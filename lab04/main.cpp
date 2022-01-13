#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

//parametry
double e = 1.0;
double delta = 0.1;
int nx=150;
int ny=100;
double V1 = 10.0; // warunek brzegowy na dole
double V2=0.0;    // warunek brzegowy na górze
double xmax=delta*nx;
double ymax=delta*ny;
double sigmaX=0.1*xmax;
double sigmaY=0.1*ymax;

double p1(double x , double y){ //gęstość
  return exp( ((-1)*((pow((x-0.35*xmax),2))/(sigmaX*sigmaX))) - (pow((y-0.5*ymax),2))/(sigmaY*sigmaY)  );
}

double p2(double x , double y){ //gęstość
  return -exp( ((-1)*((pow((x-0.65*xmax),2))/(sigmaX*sigmaX))) - (pow((y-0.5*ymax),2))/(sigmaY*sigmaY)  );
}

double S_fun(double** V, double** p){ //zastąpenie całkowania sumowaniem po węzłach
  double s=0.0;

  for(int i=0; i<=nx-1;i++ )
    for(int j=0; j<=ny-1; j++)
      s+=delta*delta*(0.5 * pow(((V[i+1][j]-V[i][j])/(delta)),2) + 0.5 * pow(((V[i][j+1]-V[i][j])/(delta)),2)- p[i][j]*V[i][j] );

  return s;
}

void RelGlob(FILE* S, FILE* V, FILE* Er, double omg){

  double **Vs, **Vn;
  Vs = new double*[nx+1];
  Vn = new double*[nx+1];
  int iter=0;
  double Sp =0.0;

  for(int i=0; i<nx+1;i++){
    Vs[i] = new double[ny+1];
    Vn[i] = new double[ny+1];
  }

  for(int i=0; i<nx+1;i++){
    for(int j=0; j<ny+1;j++)
      Vs[i][j] = Vn[i][j] = 0.0;  //załozenie zadania
  }
      
  double **p = new double*[nx];

  for(int i=0; i<nx;i++)
    p[i] = new double[ny];


  for(int i=0; i<nx;i++)
    for(int j=0; j<ny;j++)
      p[i][j] = p1(i*delta, j*delta) + p2(i*delta, j*delta); 


  for(int i=0; i<nx+1;i++){
    Vs[i][0] = Vn[i][0]=V1;   //WB dół
    Vs[i][ny] = Vn[i][ny]=V2; //WB góra
  }

  double Sc = S_fun(Vs, p);

  while(1){
    Sp = Sc;

    for(int i=1; i<nx;i++)
      for(int j=1; j<ny;j++)
        Vn[i][j] = 0.25*(Vs[i+1][j]+Vs[i-1][j]+Vs[i][j+1]+Vs[i][j-1]+(delta*delta/e)*p[i][j]);

    for(int i=0; i<ny;i++){
      Vn[0][i]=Vn[1][i];
      Vn[nx][i] = Vn[nx-1][i];
    }

    for(int i=0; i<nx+1;i++)
      for(int j=0; j<ny+1;j++)
        Vs[i][j]=(1-omg)*Vs[i][j]+ omg*Vn[i][j];

    Sc=S_fun(Vn, p);

    fprintf(S, "%d\t%g\n",iter, Sc);

    if(fabs((Sc-Sp)/Sp)<1e-8)
      break;
      
    iter++;
  }

  for(int i=0; i<nx+1;i++){
      for(int j=0; j<ny+1;j++)
        fprintf(V, "%d\t%d\t%g\n",i,j,Vn[i][j]);
      fprintf(V, "\n");
  }

  for(int i=1; i<nx;i++){
      for(int j=1; j<ny;j++){
        double err=(Vn[i+1][j]-2*Vn[i][j]+Vn[i-1][j])/(delta*delta) + (Vn[i][j+1]-2*Vn[i][j]+Vn[i][j-1])/(delta*delta) + p[i][j]/e;
        
        fprintf(Er, "%g\t%g\t%g\n",i*delta,j*delta,err);
      }
      fprintf(Er, "\n");
  }

    //zwalnianie pamieci
  for(int i=0; i<nx;i++){
    delete[] Vs[i];
    delete[] Vn[i];
    delete[] p[i];
  }
  delete[] Vn[nx]; delete[] Vs[nx]; //ostatni element
  delete[] Vn; delete[] Vs; delete[] p;

}

void RelLoc(FILE *f, double omg){

  double **V = new double*[nx+1];
  double **p = new double*[nx];
  int iter=0;
  double Sp =0.0;

  for(int i=0; i<nx+1;i++)
    V[i] = new double[ny+1];
  
  for(int i=0; i<nx+1;i++){
    for(int j=0; j<ny+1;j++)
      V[i][j] = 0.0;  //załozenie zadania
  }

  for(int i=0; i<nx;i++)
    p[i] = new double[ny];
  
  for(int i=0; i<nx;i++)
    for(int j=0; j<ny;j++)
      p[i][j] = p1(i*delta, j*delta) + p2(i*delta, j*delta); 


  for(int i=0; i<nx+1;i++){
    V[i][0]=V1;  //WB dół
    V[i][ny]=V2; //WB góra
  }

  double Sc = S_fun(V, p);

  while(1){
    Sp = Sc;

    for(int i=1; i<nx;i++)
      for(int j=1; j<ny;j++)
        V[i][j] = (1-omg) *V[i][j]+(omg/4)*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]+(delta*delta/e)*p[i][j]);

    for(int i=0; i<ny;i++){
      V[0][i]=V[1][i];
      V[nx][i] = V[nx-1][i];
    }

    Sc=S_fun(V, p);

    fprintf(f, "%d\t%g\n",iter, Sc);

    if(fabs((Sc-Sp)/Sp)<1e-8)
      break;
      
    iter++;
  }
    
}


int main(){
  FILE *Global_S_06 = fopen("Global_S_06.dat", "w");
  FILE *Global_V_06 = fopen("Global_V_06.dat", "w");
  FILE *Global_Er_06 = fopen("Global_Er_06.dat", "w");
  if(Global_S_06 && Global_V_06 && Global_Er_06) RelGlob(Global_S_06,Global_V_06,Global_Er_06, 0.6);

  fclose(Global_S_06); fclose(Global_V_06); fclose(Global_Er_06);

  FILE *Global_S_1 = fopen("Global_S_1.dat", "w");
  FILE *Global_V_1 = fopen("Global_V_1.dat", "w");
  FILE *Global_Er_1 = fopen("Global_Er_1.dat", "w");
  if(Global_S_1 && Global_V_1 && Global_Er_1) RelGlob(Global_S_1, Global_V_1, Global_Er_1, 1.0);
  
  fclose(Global_S_1); fclose(Global_V_1); fclose(Global_Er_1);


  FILE *Local_w_1 = fopen("Local_w_1.dat", "w");
  if(Local_w_1) RelLoc(Local_w_1, 1.0);
  fclose(Local_w_1);

  FILE *Local_w_14 = fopen("Local_w_14.dat", "w");
  if(Local_w_14) RelLoc(Local_w_14, 1.4);
  fclose(Local_w_14);
  
  FILE *Local_w_18 = fopen("Local_w_18.dat", "w");
  if(Local_w_18) RelLoc(Local_w_18, 1.8);
  fclose(Local_w_18);

  FILE *Local_w_19 = fopen("Local_w_19.dat", "w");
  if(Local_w_19) RelLoc(Local_w_19, 1.9);
  fclose(Local_w_19);

  return 0;
}