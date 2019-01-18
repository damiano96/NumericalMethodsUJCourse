#include <iostream>
#include <cstdio>

#define N 1000
#define h 0.01 
#define h_square h*h

void initMatrix(double *c, double *b, double *a, double *f) {
    b[0] = 1.0;
    c[0] = 0;
    f[0] = 1.0;
    for(int i=0; i<N; i++) {
        a[i] = 1.0;
        if(i>0) {
            b[i] = (h_square - 2.0);
            c[i] = 1;
            f[i] = 0;
        }
    }
    f[N] = -1.0;
    c[N-1] = 1.0;
    b[N] = -2.0;
}
void solve(double *a, double *b, double *c, double *y, double *f) {
    for(int i=1; i<=N; i++) {
        double l = a[i] / b[i-1];
        a[i] = 0;
        b[i] -= l * c[i-1];
        f[i] -= l * f[i-1];
    }
    y[N] = f[N] / b[N];
    for(int i=N-1; i>=0; i--) {
        y[i] = (f[i] - c[i] * y[i+1])/b[i];
    }
}
void printResult(double *y) {
    for(int n=0; n<=N; n++) printf("%.10lf\t%.10lf\n", n*h, y[n]);
}

int main(int argc, char const *argv[]) {
    double b[N+1], c[N], a[N], y[N+1], f[N+1];

    initMatrix(c, b, a, f);
    solve(a,b,c,y,f);
    printResult(y);

    return 0;
}
