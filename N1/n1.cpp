#include <iostream>
#include <cmath>

long double derective1(double h, double x) {
    return fabs(((cos(x+h) - cos(x)) / h) + sin(x));
}
long double derective2(double h, double x) {
    return fabs(((cos(x+h) - cos(x-h)) / (2*h)) + sin(x));
}
long double derective3(double h, double x) {
    return fabs(((-cos(x+2*h) + 8*cos(x+h) - 8*cos(x-h) + cos(x-2*h)) / (12*h)) + sin(x));
}
int main() {
    double x = 1.0;
    double h = 1.0;

    for(int i=14; i>=0; i--) {
        printf("%.22LF\t", derective1(h,x));
        printf("%.22LF \t", derective2(h,x));
        printf("%.22LF \t\t", derective3(h,x));
        printf("h = %G \n", h);
        h*=0.1;
    }
}

