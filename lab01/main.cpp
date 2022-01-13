#include <iostream>
#include <cmath>
#include <string>

double analityczne(double t, double lambda){
	return exp(t*lambda);  //y(t)
}

void Euler(double delta){

	std::string p = "Euler" + std::to_string(delta) + ".dat";
	FILE *f = fopen(&p[0], "w");

//parametry
  double y0 = 1.;
	double lambda = -1.;
  double pocz = 0.;
  double kon = 5.;

	double y = y0;
	double yp=0.;
	
	for (double i=pocz; i<=kon; i+=delta) {
		fprintf(f, "%g\t%g\t%g\t%g\n", i, y, analityczne(i, lambda), y - analityczne(i, lambda));
		yp = y;
		y = yp + delta * lambda * yp;
	}
	
	fclose(f);
}

void RK2(double delta){
  
//parametry
  double y0 = 1.;
	double lambda = -1.;
  double pocz = 0.;
  double kon = 5.;

	std::string p = "RK2" + std::to_string(delta) + ".dat";
	FILE *f = fopen(&p[0], "w");
	
	double y = y0;
	double yp=0.0;
  double k1, k2;
	
	for (double i=pocz; i<=kon; i+=delta) {
		fprintf(f, "%g\t%g\t%g\t%g\n", i, y, analityczne(i, lambda),  y - analityczne(i, lambda)) ;
		yp = y;
		k1 = lambda * yp;
		k2 = lambda * (yp +delta*k1);
		y = yp + (delta/2)* (k1 + k2);
	}
	
	fclose(f);	
}

void RK4(double delta ){

	std::string p = "RK4" + std::to_string(delta) + ".dat";
	FILE *f = fopen(&p[0], "w");
	
//parametry
  double y0 = 1.;
	double lambda = -1.;
  double pocz = 0.;
  double kon = 5.;

	double y = y0;
	double yp=0.0;
  double k1, k2, k3, k4;
	
	for (double i=pocz; i<=kon; i+=delta) {
		fprintf(f, "%g\t%g\t%g\t%g\n", i, y, analityczne(i, lambda), y-analityczne(i, lambda));
		yp = y;
		k1 = lambda * yp;
		k2 = lambda * (yp +(delta * k1/2 ) );
		k3 = lambda * (yp +(delta * k2/2 ) );
		k4 = lambda * (yp +(delta * k3) );
		y = yp + (delta/6) *(k1 +2 * k2 + 2 * k3 + k4);
	}
	
	fclose(f);	
}

//potencjal
double V( double w_V, double t){
	return 10 * sin(w_V * t);
}

void RLC(double m){
  double delta = pow( 10., -4. );
  double R = 100, L = 0.1, C = 0.001;
  
  double Q_0 = 0., I_0 = 0.;
  
	double w_0 = 1 / sqrt( L * C );
	double T_0 = 2 * M_PI / w_0;
	double w_V = m * w_0;
	double k1_Q, k2_Q, k3_Q, k4_Q, k1_I, k2_I, k3_I, k4_I;
	double Q = Q_0, I = I_0;
	
	std::string filename = "RLC" + std::to_string(m) + ".dat";
	FILE *file = fopen(&filename[0], "w");
	
	
	for (double t = 0.0; t <= 4 * T_0 ; t+=delta){
		k1_Q = I;
		k1_I = V( w_V, t)/L - ( 1 / ( L * C ) ) * Q - ( R / L ) * I;
		k2_Q = I + ( delta / 2 ) * k1_I;
		k2_I = V( w_V, t + (delta / 2) ) / L - ( 1 / ( L * C ) ) * ( Q + ( delta / 2 ) * k1_Q ) - ( R / L ) * ( I + ( delta / 2 ) * k1_I );
		k3_Q = I + ( delta / 2 ) * k2_I;
		k3_I = V( w_V, t + (delta / 2) ) / L - ( 1 / ( L * C ) ) * ( Q + ( delta / 2 ) * k2_Q ) - ( R / L ) * ( I + ( delta / 2 ) * k2_I );
		k4_Q = I + delta * k3_I;
		k4_I = V( w_V, t + delta ) / L - ( 1 / ( L * C ) ) * ( Q + delta * k3_Q ) - ( R / L ) * ( I + delta * k3_I ); 
	
		Q += ( delta / 6 )*( k1_Q + 2* k2_Q + 2 * k3_Q + k4_Q );
		I += ( delta / 6 ) *( k1_I + 2* k2_I + 2 * k3_I + k4_I );

		fprintf( file, "%g\t%g\t%g\n", t, I, Q );
	}
	
	fclose(file);	
}


int main(){

  double delta_t1 = 0.01;
  double delta_t2 = 0.1;
  double delta_t3 = 1.0;

	Euler(delta_t1);
	Euler(delta_t2);
	Euler(delta_t3);
	
	RK2(delta_t1);
	RK2(delta_t2);
	RK2(delta_t3);

	RK4(delta_t1);
	RK4(delta_t2);
	RK4(delta_t3);

  double w1 = 0.5; RLC(w1);
  double w2 = 0.8; RLC(w2);
  double w3 = 1.0; RLC(w3);
  double w4 = 1.2; RLC(w4);

  return 0;
}
