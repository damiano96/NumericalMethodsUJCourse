#include <iostream>
#include <cmath>
#include <iomanip>

#define E 1e-6

double f(double x) {
	return exp(cos(x) * cos(x));
}
double method_Trapez(double a, double b, double (*f)(double)) {
    double sum = f(a) + f(b);
    double prevIntegral, integral = 0.0; 
    int intervals = 1;
    do  {
        intervals++;
        prevIntegral = integral;
        double sum2 = 0.0;
 	    double h = (b-a) / intervals;

        for(int i=1; i<intervals; i++) sum2 += f(a+i*h);

	    integral = h*((0.5 * sum) + sum2);
    } while (fabs(integral - prevIntegral) > E);
    std::cout << "Metoda trapezow - liczba iteracji: " << intervals << " - \t";
    return integral;
}
double method_Simpson(double a, double b, double (*f)(double)) {
    double sum = f(a) + f(b);
    double prevIntegral, integral = 0.0; 
    int intervals = 0;
    do {
        intervals++;
        prevIntegral = integral;
        double sum2 = 0.0;
        double sum4 = 0.0;
        double h = (b-a)/intervals;
        double h2 = h/2.0;

        for(int i=1; i<=intervals; i++) sum4 += f(a+i*h-h2);
        for(int i=1; i<intervals; i++) sum2 += f(a+i*h);
        
        integral = h/6.0 * (sum + 4.0 * sum4 + 2.0 * sum2 );
    } while (fabs(integral - prevIntegral) > E);
    std::cout << "Metoda Simpsona - liczba iteracji: " << intervals << " - \t";
	return integral;
}
double method_threeAndEight(double a, double b, double (*f)(double)) { 
	double sum = f(a) + f(b);
    double sum_copy = sum;
    double prevIntegral, integral = 0.0; 
    int intervals = 0;
    do {
        intervals+=3;
        sum = sum_copy;
        prevIntegral = integral;
        double h = (b-a) / intervals;
        for (int i=1; i<intervals; i++) {
            if(i%3 == 0) sum += 2*f(a+i*h);
            else sum += 3*f(a+i*h);
        }
        integral = (3.0*h / 8.0) * sum;
        
    } while (fabs(integral - prevIntegral) > E);
    std::cout << "Metoda 3/8 - liczba iteracji: " << intervals << " - \t";
	return integral;
}

int main() {
    std::cout << std::fixed << std::setprecision(10) << method_Trapez(0, M_PI, f) << std::endl;
    std::cout << std::fixed << std::setprecision(10) << method_Simpson(0, M_PI, f) << std::endl;
    std::cout << std::fixed << std::setprecision(10) << method_threeAndEight(0, M_PI, f) << std::endl;
}