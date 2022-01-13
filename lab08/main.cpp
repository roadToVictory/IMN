#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

//parametry
#define delta 0.01
#define nx 400
#define ny 90
#define i1 200
#define i2 210
#define j1 50
#define sigma 10*delta
#define xa .45
#define ya .45

void solver(double D){
    std::ifstream read; read.open("psi.dat");
    if(!read)
        return;
    
    double** Psi = new double*[nx+1];
    double** Vx = new double*[nx+1];
    double** Vy = new double*[nx+1];
    double** u0 = new double*[nx+1];
    double** u1 = new double*[nx+1];
    double deltaT=.0;

    for(int i=0; i<=nx; i++){
        Psi[i] = new double[ny+1];
        Vx[i] = new double[ny+1];
        Vy[i] = new double[ny+1];
        u0[i] = new double[ny+1];
        u1[i] = new double[ny+1];
    }
        

    for(int i=0; i<=nx; i++)
        for(int j=0; j<=ny; j++)
            Psi[i][j] = Vx[i][j] = Vy[i][j] = u0[i][j] = u1[i][j] =.0;


    int i,j;
    double x;
    while(read>>i>>j>>x)
        Psi[i][j] = x;
    read.close();

    //pole prędkosci 
    for(int i=1; i<=nx-1; i++)
        for(int j=1; j<=ny-1; j++){
            Vx[i][j] = (Psi[i][j+1]-Psi[i][j-1])/(2*delta);
            Vy[i][j] = -(Psi[i+1][j]-Psi[i-1][j])/(2*delta);
        }
    
    for(int i=i1; i<=i2; i++)          //na zastawce
        for(int j=0; j<=j1; j++)
            Vx[i][j] = Vy[i][j] = .0;
        
    for(int i=1; i<=nx-1; i++)		   //na górnym i dolnym brzegu
       Vx[i][0] = Vy[i][ny] = 0.0;
    
        
    for(int j=0; j<=ny; j++) {          //na lewyn i prawym
        Vx[0][j] = Vx[1][j];
        Vx[nx][j] = Vx[nx-1][j];
    }

    FILE *f1= fopen(("V"+std::to_string(D)+".dat").c_str(), "w");
    if(!f1) return;

    for(int i=0; i<=nx; i++) {   
      for(int j=0; j<=ny; j++)
          fprintf(f1, "%g\t%g\t%g\t%g\n", i*delta, j*delta, Vx[i][j], Vy[i][j]);
        fprintf(f1, "\n");
    } fclose(f1);

    double Vmax=sqrt(pow(Vx[0][0],2)+ pow(Vy[0][0],2));
    double tmp=.0;
    for(int i=0; i<=nx; i++) 
      for(int j=0; j<=ny; j++){
          tmp = sqrt(pow(Vx[i][j],2)+ pow(Vy[i][j],2));
          Vmax = (tmp>Vmax)? tmp : Vmax ; 
      }

    deltaT = delta/(4*Vmax);    //krok czasowy

            //algorytm dla rownania AD
    for(int i=0; i<=nx; i++)
        for(int j=0;  j<=ny; j++)
            u0[i][j] = 1./(2*M_PI*pow(sigma,2))*exp(-( pow((delta*i-xa),2) + pow((delta*j-ya),2)  )/(2*pow(sigma,2)));

    for(int it=1; it<=10000; it++){
        for(int i=0; i<=nx; i++)
            for(int j=0;  j<=ny; j++)
                u1[i][j] = u0[i][j];

        for(int k=1; k<=20; k++){
            for(int i=0; i<=nx; i++){
                for(int j=1;  j<=ny-1; j++){
                    if(i>=i1 && i<=i2 && j<=j1){
                        continue;}
                    else if(i==0){
                        u1[i][j] = (1./(1.+(2*D*deltaT)/(pow(delta,2))))*( u0[i][j] - (deltaT/2 *Vx[i][j]) * ( (u0[i+1][j]-u0[nx][j])/(2*delta) + (u1[i+1][j]-u1[nx][j])/(2*delta) ) 
                        - (deltaT/2 * Vy[i][j]) * ((u0[i][j+1]-u0[i][j-1])/(2*delta) + (u1[i][j+1]-u1[i][j-1])/(2*delta)) 
                        +  (deltaT/2)*D* (  (u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j])/(pow(delta,2))  + ( (u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1])/(pow(delta,2))   )  ) );
                    }
                    else if(i==nx){
                        u1[i][j] = (1./(1.+(2*D*deltaT)/(pow(delta,2))))*( u0[i][j] - (deltaT/2 *Vx[i][j]) * ( (u0[0][j]-u0[i-1][j])/(2*delta) + (u1[0][j]-u1[i-1][j])/(2*delta) ) 
                        - (deltaT/2 * Vy[i][j]) * ((u0[i][j+1]-u0[i][j-1])/(2*delta) + (u1[i][j+1]-u1[i][j-1])/(2*delta)) 
                        +  (deltaT/2)*D* (  (u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j])/(pow(delta,2))  + ( (u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1])/(pow(delta,2))   )  ) );
                    }
                    else{
                         u1[i][j] = (1./(1.+(2*D*deltaT)/(pow(delta,2))))*( u0[i][j] - (deltaT/2 *Vx[i][j]) * ( (u0[i+1][j]-u0[i-1][j])/(2*delta) + (u1[i+1][j]-u1[i-1][j])/(2*delta) ) 
                        - (deltaT/2 * Vy[i][j]) * ((u0[i][j+1]-u0[i][j-1])/(2*delta) + (u1[i][j+1]-u1[i][j-1])/(2*delta)) 
                        +  (deltaT/2)*D* (  (u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j])/(pow(delta,2))  + ( (u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1])/(pow(delta,2))   )  ) );
                    }
                }
            }
        }

        for(int i=0; i<=nx; i++)
            for(int j=0;  j<=ny; j++)
                u0[i][j] = u1[i][j];

        double c=.0, xsr=.0;
        for(int i=0; i<=nx; i++)
            for(int j=1;  j<=ny; j++){
                c+=u0[i][j];
                xsr+=i*delta*u0[i][j];
            }
        c*=pow(delta,2);
        xsr*=pow(delta,2);

        FILE *f2= fopen(("C"+std::to_string(D)+".dat").c_str(), "w");
        if(!f2) return;
        FILE *f3= fopen(("x"+std::to_string(D)+".dat").c_str(), "w");
        if(!f3) return;
        fprintf(f2, "%g\t%g\n", delta*it, c);
        fprintf(f3, "%g\t%g\n", delta*it, xsr);
        fclose(f2); fclose(f3);
        
        FILE *f4= fopen(("u"+std::to_string(D)+".dat").c_str(), "w");
        if(!f4) return;

            for(int i=0; i<=nx; i++){
                for(int j=0;  j<=ny; j++)
                    fprintf(f4, "%d\t%d\t%g\n", i, j, u1[i][j]);
                fprintf(f4, "\n");
            }
        fclose(f4);
    }

    for(int i=0; i<=nx; i++){
        delete[] Psi[i];
        delete[] Vx[i];
        delete[] Vy[i];
        delete[] u0[i];
        delete[] u1[i];
    }
    delete[] Psi; delete[] Vx; delete[] Vy; delete[] u0; delete[] u1;
}

int main(){
  solver(0.);
 // solver(.1);
}