#include <iostream>
#include <cmath>
#include <string>

//parametry
double beta = 0.001;
int N = 500, tmax=100;
double gama = 0.1;
double delta=0.1;
double u0=1;
double TOL=10e-6;
double alfa= beta*N-gama;

void Picard(){
  FILE *p = fopen("picard.dat", "w");
  if(p){
    double un = u0;
    double un1 = 0.0;
    double unm = u0;

    for(double i=0.; i<=tmax; i+=delta){
      for(int mi=0; mi<=20; mi++){
        un1=un+(delta/2)*((alfa*un-beta*un*un)+(alfa*unm-beta*unm*unm));
        if(fabs(un1-unm) <TOL) break;
        unm=un1; //poprzednie = kolejne -> do nastepnej iteracji
      }
      un=unm=un1;
      fprintf(p, "%g\t%g\t%g\n", i, un, N-un);
    }
    fclose(p);
  }
  else perror("file open error");
}

void Newton(){
  FILE *p = fopen("Newton.dat", "w");
  if(p){
    double un = u0;
    double un1 = 0.0;
    double unm = u0;

    for(double i=0.; i<=tmax; i+=delta){
      for(int mi=0; mi<=20; mi++){
        un1= unm-(unm-un-(delta/2)*(alfa*un-beta*un*un+(alfa*unm-beta*unm*unm)))/(1-(delta/2)*(alfa-2*beta*unm));
        if(fabs(un1-unm) <TOL) break;
        unm=un1; //poprzednie = kolejne -> do nastepnej iteracji
      }
      un=unm=un1;
      fprintf(p, "%g\t%g\t%g\n", i, un, N-un);
    }

    fclose(p);
  }
  else perror("file open error");
}

double f(double u){
  return (beta*N-gama)*u-beta*u*u;
}

void RK2(){
  FILE *p = fopen("RK2.dat", "w");
  if(p){
    double a[2][2]={0.25, 0.25-sqrt(3)/6, 0.25+sqrt(3)/6, 0.25};
    double b1=.5, b2=.5, U1=.0, U2=.0, dU1=0., dU2=0., m[2][2], F1, F2, un=u0, un1=.0, U1m, U2m ;

    for(double i=.0; i<=tmax; i+=delta){
      U1=U2=un;
      m[0][0]=1-delta*a[0][0]*(alfa-2*beta*U1);
      m[0][1]=-delta*a[0][1]*(alfa-2*beta*U2);
      m[1][0] = -delta*a[1][0]*(alfa-2*beta*U1);
      m[1][1]=1-delta*a[1][1]*(alfa-2*beta*U2);

      F1 = U1-un-delta*(a[0][0]*(alfa*U1-beta*U1*U1)+a[0][1]*(alfa*U2-beta*U2*U2));
      F2 = U2-un-delta*(a[1][0]*(alfa*U1-beta*U1*U1)+a[1][1]*(alfa*U2-beta*U2*U2));

      dU1= (F2*m[0][1]-F1*m[1][1])/(m[0][0]*m[1][1]-m[0][1]*m[1][0]);
      dU1= (F1*m[1][0]-F2*m[0][0])/(m[0][0]*m[1][1]-m[0][1]*m[1][0]);

      for(int mi=0; mi<=20; mi++){
        U1m = U1+dU1;
        U2m=U2+dU2;
        if(fabs(U1m-U1)<TOL || fabs(U2m-U2)<TOL) break;
        else{ U1=U1m; U2=U2m;}
      }

      un1 = un+delta*(b1*f(U1)+b2*f(U2));
      un=U1=U2=un1;
      fprintf(p, "%g\t%g\t%g\n", i, un, N-un);
    }
    fclose(p);
  }
  else perror("file open error");
}

int main(){
//metoda trapezów z iteracją Picarda
  Picard();
//metoda trapezów z iteracją Newtona
  Newton();
//niejawna RK2
  RK2();
}