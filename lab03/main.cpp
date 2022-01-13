#include <iostream>
#include <cmath>
#include <string>
#include <utility>
#include <functional>

//parametry zadania
double x0=0.01;
double v0=0.;
double delta_t0=1.0;
double S=0.75;
double p=2;
double tmax=40.;
double alfa=5.0;

using pairOfDouble = std::pair<double,double>; 

double g(double x ,double v){return alfa*(1-x*x)*v-x;}

double max(double a, double b){return a>b?a:b;}

pairOfDouble trapez(double xn, double vn, double dt, double alpha){
  
  double xn1=xn, vn1=vn, s=1e-10, dx=0., dv=0.;
  do{
    double a[2][2] = {1.0, -dt/2, (-dt/2)*(-2*alpha*xn1*vn1-1) , 1-(dt/2)*alpha*(1-pow(xn1,2)) };

    double F = xn1-xn-(dt/2)*(vn+vn1);
    double G = vn1-vn-(dt/2)*(g(xn,vn)+g(xn1,vn1));
    dx=((-F)*a[1][1]-(-G)*a[0][1])/(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    dv=(a[0][0]*(-G)-a[1][0]*(-F))/(a[0][0]*a[1][1]-a[0][1]*a[1][0]);

    xn1+=dx;
    vn1+=dv;
  }while(fabs(dx)<s && fabs(dv)<s);

  return pairOfDouble(xn1,vn1);
}

pairOfDouble RK2(double xn, double vn, double dt, double alpha){

  double k1x = vn;
  double k1v = g(xn,vn);
  double k2x=vn+dt*k1v;
  double k2v=alpha*(1-pow(xn+dt*k1x,2))*(vn+dt*k1v)-(xn+dt*k1x);
  double xn1 = xn+(dt/2)*(k1x+k2x);
  double vn1= vn+(dt/2)*(k1v+k2v);

  return pairOfDouble(xn1,vn1);
}

void step(FILE* f, pairOfDouble(*fun)(double xn, double vn, double dt, double alpha), double TOL){
  double t=0., dt=delta_t0, xn=x0, vn=v0;

  pairOfDouble single(xn,vn); //krok pojedynczy
  pairOfDouble dual(xn,vn); //krok podwójny
  double Ex=0.0, Ev=0.0;

  do{
    single = fun(xn,vn,dt,alfa);
    single = fun(single.first, single.second, dt, alfa);

    dual=fun(xn,vn, 2*dt, alfa);

    Ex= (single.first-dual.first)/(pow(2,p)-1);
    Ev= (single.second-dual.second)/(pow(2,p)-1);

    if(max(fabs(Ex), fabs(Ev))<TOL){
      t+=2*dt;
      xn=single.first;
      vn=single.second;
      fprintf(f, "%g\t%g\t%g\t%g\n", t, dt, xn,vn);
    }

    dt*=pow(((S*TOL)/(max(fabs(Ex), fabs(Ev)))), (1./(p+1.)));

  }while(t<tmax);
}


int main(){

  double TOL2 = 1e-2;
  double TOL5 = 1e-5;
  
  FILE* trapezTOL2 = fopen("tTOL2.dat", "w");
  FILE* trapezTOL5 = fopen("tTOL5.dat", "w");
  FILE* rk2TOL2 = fopen("rTOL2.dat", "w");
  FILE* rk2TOL5 = fopen("rTOL5.dat", "w");
  
  if(trapezTOL2) step(trapezTOL2, &trapez, TOL2);
  if(trapezTOL5) step(trapezTOL5, &trapez, TOL5);

  if(rk2TOL2) step(rk2TOL2, &RK2, TOL2);
  if(rk2TOL5) step(rk2TOL5, &RK2, TOL5);

  fclose(trapezTOL2); fclose(trapezTOL5); fclose(rk2TOL2); fclose(rk2TOL5);
}