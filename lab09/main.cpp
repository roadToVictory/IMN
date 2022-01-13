#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

#include</usr/include/gsl/gsl_math.h>
#include</usr/include/gsl/gsl_linalg.h>
#include</usr/include/gsl/gsl_blas.h>

#define nx 40
#define ny 40
#define N (nx+1)*(ny+1)
#define delta 1
#define deltaT 1
#define Ta 40
#define Tb 0
#define Tc 30
#define Td 0
#define kb 0.1
#define kD 0.6
#define IT_MAX 2000

int main(){
    gsl_matrix* A = gsl_matrix_calloc(N,N);
    gsl_matrix* B = gsl_matrix_calloc(N,N);
    gsl_vector* c = gsl_vector_calloc(N);
    gsl_vector* T = gsl_vector_calloc(N);
    gsl_vector* d = gsl_vector_calloc(N);
    gsl_permutation* p = gsl_permutation_alloc(N);
    int signum=.0;

    int l=0; double val=.0;

    for(int i=1; i<=nx-1;i++){  //wnetrze
        for(int j=1; j<=ny-1;j++){
            l=i+j*(nx+1);
            val = deltaT/(2*pow(delta,2));
            gsl_matrix_set(A, l, l-nx-1, val );
            gsl_matrix_set(A, l, l-1, val );
            gsl_matrix_set(A, l, l+1, val );
            gsl_matrix_set(A, l, l+nx+1, val );
            gsl_matrix_set(A, l, l, (-4)*val-1 );
            gsl_matrix_set(B, l, l-nx-1, -val );
            gsl_matrix_set(B, l, l-1, -val );
            gsl_matrix_set(B, l, l+1, -val );
            gsl_matrix_set(B, l, l+nx+1, -val );
            gsl_matrix_set(B, l, l, 4*val-1 );
        }
    }

    for(int i=0; i<=nx;i+=nx){  //warunki Dirichleta -> lewy i prawy
        for(int j=0; j<=ny;j++){
            l=i+j*(nx+1);   //i=0, i=nx
            gsl_matrix_set(A, l, l, 1 );
            gsl_matrix_set(B, l, l, 1 );
            gsl_vector_set(c, l, 0);
        }
    }

     for(int i=1; i<=nx-1;i++){  //warunki von Neumanna -> górny dla n+1
        l=i+ny*(nx+1);  
        gsl_matrix_set(A, l, l-nx-1, -1/(kb*delta) );
        gsl_matrix_set(A, l, l, 1+(1/(kb*delta)) );
        gsl_vector_set(c, l, Tb);

        for(int j=0; j<=ny; j++)
            gsl_matrix_set(B, l, j, 0 );
    }

    for(int i=1; i<=nx-1;i++){  //warunki von Neumanna -> dolny dla n+1
        l=i; 
        gsl_matrix_set(A, l, l+nx+1, -1/(kD*delta) );
        gsl_matrix_set(A, l, l, 1+(1/(kD*delta)) );
        gsl_vector_set(c, l, Td);

        for(int j=0; j<=ny; j++)
            gsl_matrix_set(B, l, j, 0 );
    }

    // for(int i=1; i<=nx;i++){ 
    //     for(int j=0; j<=ny;j++){
    //         l=i+j*(nx+1);
    //         gsl_vector_set(T,l,0);
    //     }
    // }

    for(int j=0; j<=ny;j++){ 
        l=j*(nx+1);
        gsl_vector_set(T,l,Ta);
    }

    for(int j=0; j<=ny;j++){ 
        l=nx+j*(nx+1);
        gsl_vector_set(T,l,Tc);
    }

    for(int i = 1; i <= nx - 1; i++){	 	// pozostały obszar
        for(int j = 0; j <= ny; j++){
			l = i + j * ( nx + 1 );
			gsl_vector_set(T, l, 0);
		}
    }

     gsl_linalg_LU_decomp(A , p , &signum );

    FILE *Temp = fopen("T.dat", "w");
	FILE *op = fopen("op.dat", "w");

    for(int it=0; it<=IT_MAX; it++){
         gsl_blas_dgemv(CblasNoTrans, 1 ,B , T ,0 , d);
         gsl_blas_daxpy(1, c, d );
         gsl_linalg_LU_solve (A ,p , d , T);

         int iter[5]={100,200,500,1000,2000};
         for(int k=0; k<5; k++){
             if(iter[k]==it){
                for(int i = 1; i <= nx-1; i++){
				    for(int j = 1; j <= ny-1; j++){
                        l=i+j*(nx+1);
                        fprintf(op, "%d\t%d\t%d\t%g\n", it, i*delta, j *delta, (((gsl_vector_get(T, l+1) - 2*gsl_vector_get(T,l) + gsl_vector_get(T, l-1))/pow(delta,2)) + ((gsl_vector_get(T, l+nx+1) - 2*gsl_vector_get(T,l) + gsl_vector_get(T, l-nx-1))/pow(delta,2))));
                        fprintf(Temp, "%d\t%d\t%d\t%g\n", it, i*delta, j *delta, gsl_vector_get(T, l));
                    }
                    fprintf(op, "\n");
                    fprintf(Temp, "\n");
                }
                fprintf(op, "\n");
                fprintf(Temp, "\n");
            }
         }
         

     }

    gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(c);
	gsl_vector_free(T);
    gsl_vector_free(d);
    gsl_permutation_free(p);
}
