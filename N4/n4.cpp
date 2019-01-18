#include <iostream>
#include <cstdio>

#define N 1000
#define h 0.01
#define h_square h*h


int main() {
    double y[N];

    y[0] = 1;
    y[1] = -4.0 / (h_square - 6.0);
    y[2] = 3.0 + (16.0/ (h_square - 6.0));

    for(int i=2; i<N-1; i++) {
        y[i+1] = y[i]*(-h_square + 2.0) - y[i-1];
    }

    for(int i=0; i<N; i++) {
        printf("%.11lf\t %.11lf\n", i*h, y[i]);
    }
}