#include <iostream>
#include <iomanip>
#include <cmath>

#define E 1e-10
int i1 = 0;
int i2 = 0;
int i3 = 0;

double f(double l) {
    return -(l*l*l) + (12 * (l * l)) - (46 * l) + 56;
}

double df(double l) {
    return -3*(l*l) + (24*l) - 46;
}

double secantMethod(double (*f)(double), double a, double b) {
	double x0 = 0.0;

	while(true) {
		double f_a = f(a);
		double f_b = f(b);
		
		if(std::fabs(f_b - f_a) < E) break;
		
		x0 = b - f_b * (b - a) / (f_b - f_a);
		a = b;
		b = x0;
        i1++;
	}
	return x0;
}

double methodNewton(double (*f)(double), double start) {
	double x0 = start;
	double f_x0 = f(x0);

	while(std::fabs(f_x0) > E) {
		double d_f = df(x0);
		x0 -= f_x0 / d_f;
		f_x0 = f(x0);
        i2++;
	}
	return x0;
}

double regulaFalsi(double (*f)(double), double a, double b) {
    double fa = f(a), fb = f(b), x1 = a, x0 = b, f0;
    while(fabs(x1 - x0) > E) {
        x1 = x0;
        x0 = a - fa * (b - a) / (fb - fa);
        f0 = f(x0);
        if(fabs(f0) < E) break;
        if(fa * f0 < 0){
            b = x0; 
            fb = f0;
        }
        else {
            a = x0; 
            fa = f0;
        }
        i3++;
    }
    return x0;
}
int main(int argc, char const *argv[]) {
    
    double x0, x1;

    printf("Regula Falsi\n");
	for(double i=2; i<6; i+=0.01) {
		x0 = i;
		x1 = 0.01 + x0;
		if(f(x0) * f(x1) <= 0) {
			std::cout << std::fixed << std::setprecision(10) << regulaFalsi(f, x0, x1) << std::endl;
		}
	} 
    std::cout << "Liczba iteracji: " << i3 << std::endl  << std::endl;

    printf("Metoda Newtona\n");
	for(double i=2; i<6; i+=0.01) {
		x0 = i;
		x1 = x0 + 0.01;
		if(f(x0) * f(x1) <= 0) {
            if(methodNewton(f, x0) >= x0 && methodNewton(f, x0) <= x1) {
			    std::cout << std::fixed << std::setprecision(10) << methodNewton(f, x0) << std::endl;
            }
		}
	} 
    std::cout << "Liczba iteracji: " << i2;

    std::cout << "\n\n";
    printf("Metoda siecznych:\n");
	for (double i=2; i<6; i+=0.01) {
		x0 = i;
		x1 = x0 + 0.01;
		if(f(x0) * f(x1) <= 0 ) {
            std::cout << std::fixed << std::setprecision(10) << secantMethod(f, x0, x1) << std::endl;
        }
	}

    std::cout << "Liczba iteracji: " << i1 << std::endl;
    
}
