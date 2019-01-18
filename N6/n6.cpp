#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#define precision 1e-10

double getActualNorm(std::vector<double> v) {
	double result = 0;
	for(int i=0; i<v.size(); i++) {
		result += v[i] * v[i];
	}
	return std::sqrt(result);
}
std::vector<double> multiplyVector(std::vector<double> A1, std::vector<double> A2, std::vector<double> A3, std::vector<double> y) {
    std::vector<double> returnVector(y.size(), 0.0);
    returnVector[0] = ((A2[0] * y[0]) + (A1[0] * y[1])); 
    returnVector[1] = ((A3[0] * y[0]) + (A2[1] * y[1]) + (A1[1] * y[2])); 
    returnVector[2] = ((A3[1] * y[1]) + (A2[2] * y[2])); 
    return returnVector;
}
double multiplyVectorT(std::vector<double> x, std::vector<double> y) {
    double result = 0;
    for(int i=0; i<x.size(); i++) {
        result += x[i] * y[i];
    }
    return result;
}
std::vector<double> multiplyVectorScalar(double x, std::vector<double> y) {
    std::vector<double> returnVector(y.size(), 0.0);
    for(int i=0; i<y.size(); i++) {
        returnVector[i] = x * y[i];
    }
    return returnVector;
}
std::vector<double> addVector(std::vector<double> x, std::vector<double> y) {
    std::vector<double> returnVector(x.size(), 0.0);
    for(int i=0; i<x.size(); i++) {
        returnVector[i] = x[i] + y[i];
    }
    return returnVector;
}
std::vector<double> substractVector(std::vector<double> x, std::vector<double> y) {
    std::vector<double> returnVector(x.size(), 0.0);
    for(int i=0; i<x.size(); i++) {
        returnVector[i] = x[i] - y[i];
    }
    return returnVector;
}

std::vector<double> solve(std::vector<double> A1, std::vector<double> A2, std::vector<double> A3, std::vector<double> b) {
    std::vector<double> r(b.size(), 0.0);
    std::vector<double> x(b.size(), 1.0);
    
    r = substractVector(b, multiplyVector(A1, A2, A3, x));
    std::vector<double> p(r);
    int iterator = 0;

    while(getActualNorm(r) > precision) {
        double rt_r = multiplyVectorT(r,r);
        std::vector<double> Ap = multiplyVector(A1, A2, A3, p);
        double alfa = rt_r / multiplyVectorT(p, Ap);

        r = substractVector(r, multiplyVectorScalar(alfa, Ap));
        double beta = multiplyVectorT(r,r) / rt_r;
        x = addVector(x, multiplyVectorScalar(alfa, p));
        p = addVector(r, multiplyVectorScalar(beta, p)); 
        iterator++;
    }
    return x;
}

int main() {
    std::vector<double> A1 {-1.0, -1.0};
    std::vector<double> A2 {4.0, 4.0, 4.0};
    std::vector<double> A3 {-1.0, -1.0};

    std::vector<double> b {2.0, 6.0, 2.0};
    std::vector<double> wynik = solve(A1, A2, A3, b);

    for(int i=0; i<wynik.size(); i++) {
        std::cout <<  std::fixed << std::setprecision(10) << wynik[i] << " \n";
    }
}