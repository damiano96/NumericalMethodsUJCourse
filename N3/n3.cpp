#include <iostream>
#include <cstdio>


#define N 1000
#define h 0.01 
#define h_square h*h

void initMatrix(double *a, double *b, double *c, double *f) {
    b[0] = 1.0;
    c[0] = 0;
    f[0] = 1.0;
    for(int i=0; i<N; i++) {
        a[i] = -1.0;
        if(i>0) {
            b[i] = (h_square + 2.0);
            c[i] = -1;
            f[i] = 0;
        }
    }
    c[N-1] = -1.0;
    b[N] = -2.0;
}
void metodaThomasa(double *a, double *b, double *c, double *u, double *f) {
	for(int i = 1; i < N; i++){
        double m = a[i] / b[i-1];
        a[i] = 0.0;
        b[i] -= m*c[i-1];
        f[i] -= m*f[i-1];
    }
    u[N] = f[N]/b[N];

    for (int k = N; k >= 0 ; k--) {
        u[k] = (f[k] - c[k]*u[k+1])/b[k];
    }
}

void solve(double *wynik, double *u, double *z, double *v) {
	double vUp = 0;
	double zDown = 0;

	for(int i = 0; i<N+1; i++){
    	vUp += v[i]*u[i];
    	zDown += v[i]*z[i];
  	}
    zDown += 1.0;
  	for(int i=0; i<N+1; i++) {
          wynik[i] = u[i] - (vUp/zDown) * z[i];
  	}
}

void printArray(double *array) {
	for (int i = 0; i <= N; i++) {
         printf("%.10f %.10f \n",  (i)*h, array[i]);
    }
}

int main() {
    double a[N], b[N+1], c[N], f[N+1], u[N+1], z[N+1], wynik[N+1];
    double u_V[N+1] = {0};
    double V[N] = {0};
    u_V[N] = 1.0;
    V[0] = 1.0;

    initMatrix(a,b,c,f);

    metodaThomasa(a,b,c,u,f);
    metodaThomasa(a,b,c,z,u_V);

    solve(wynik, u, z, V);

    printArray(wynik);

    
}