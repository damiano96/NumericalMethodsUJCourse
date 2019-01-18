#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

constexpr unsigned N = 1001;
constexpr double h = 0.01;
constexpr double h_sqare = h*h;

double errorNorm(std::vector<double>error) {
    double norm = 0.0;
    for(int i=0; i<error.size(); i++) {
        norm += error[i]*error[i];
    }
    return std::sqrt(norm);
}

std::vector<double> errorResult(std::vector<double> a, std::vector<double> b) {
    std::vector<double> result(N);
    for(int i=0; i<N; i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

void initMatrix(double *lower, double *main, double *upper, double *b) {
    main[0] = b[0] = 1.0;
    upper[0] = 0;
    lower[0] = -1.0;
    for(int i=1; i<N; i++) {
        main[i] = h_sqare + 2;
        upper[i] = -1.0;
        b[i] = 0;
        if(i<N-2) lower[i] = -1.0;
    }
    b[N-1] = -1.0;
    lower[N-2] = 1.0;
    main[N-1] = -2.0;
}

std::vector<double> jacobi_Method(double *lower, double *main, double *upper, double *b, double precision) {
    std::vector<double> x(N, 0.0);
    std::vector<double> x_new(N, 0.0);
    std::vector<double> error(N, 0.0);
    double actualPrecisionValue = 1.0;
    int iteration = 0;

    while(actualPrecisionValue >= precision) {
        x_new[0] = b[0];
        for(int i=1; i<N-1; i++)   x_new[i] = (1/main[i]) * (b[i] - (lower[i-1] * x[i-1] + upper[i] * x[i+1]));
        
        x_new[N-1] = (1/main[N-1]) * (b[N-1] - (lower[N-2] * x[N-2]));

        error = errorResult(x, x_new);
        actualPrecisionValue = errorNorm(error);
        x = x_new;
        iteration++;
    } 
    //std::cout << iteration << std::endl;

    return x;
}
std::vector<double> gauss_Seidel_Method(double *lower, double *main, double *upper, double *b, double precision) {
    std::vector<double> x(N, 0.0);
    std::vector<double> x_new(N, 0.0);
    std::vector<double> error(N, 0.0);
    double actualPrecisionValue = 1.0;
    int iteration = 0;

    while(actualPrecisionValue >= precision) {
        x_new[0] = b[0];
        x_new[N-1] =  (1/main[N-1]) * (b[N-1] - lower[N-2] * x_new[N-2]);

        for(int i=1; i<N-1; i++)    x_new[i] = (1/main[i]) * (b[i] - lower[i-1] * x_new[i-1] - upper[i] * x[i+1]);
        
        error = errorResult(x, x_new);
        actualPrecisionValue = errorNorm(error);
        x = x_new;
        iteration++;
    }
    //std::cout << iteration << std::endl;
    return x;
}
std::vector<double> successive_OverRelaxation_Method (double *lower, double *main, double *upper, double *b, double precision) {
    double w_optimal = 2.0 / (1.0 + sin(M_PI * (1.0 / (N + 1))));
    std::vector<double> x(N, 0.0);
    std::vector<double> x_new(N, 0.0);
    std::vector<double> error(N, 0.0);
    double actualPrecisionValue = 1.0;
    int iteration = 0;
    while(actualPrecisionValue >= precision) {
        x_new[0] = ((1.0 - w_optimal) * x[0] + (w_optimal * b[0]));
        for(int i=1; i<N-1; i++) {
            double x_new_Gauss = ((1.0 / main[i]) * (b[i] - lower[i-1] * x_new[i-1] - upper[i] * x[i+1])); 
            x_new[i] = ((1.0 - w_optimal) * x[i]) + w_optimal * x_new_Gauss;
        }
        x_new[N-1] = ((1.0 - w_optimal) * x[N-1]) + ((w_optimal / main[N-1]) * (b[N-1] - lower[N-2] * x_new[N-2]));

        error = errorResult(x, x_new);
        actualPrecisionValue = errorNorm(error);
        x = x_new;
        iteration++;
    } 
    //std::cout << iteration << std::endl;
    return x;
}
std::vector<double>relaxation_Richardson_Method(double A[][3], double *b, double precision) {
    std::vector<double> x(3, 1.0);
    std::vector<double> x_new(3, 0.0);
    std::vector<double> tempValues(3, 0.0);
    std::vector<double> error(3, 0.0);    
    double omega = 0.25;
    double actualPrecisionValue = 1.0;
    int iteration = 0;

    while(actualPrecisionValue >= precision) {
        for(int i=0; i<3; i++)     tempValues[i] = b[i] - ((A[i][0] * x[0]) + (A[i][1] * x[1]) + (A[i][2] * x[2]));
        
        for(int i=0; i<3; i++)     x_new[i] = x[i] + (omega * tempValues[i]);
        
        for(int i=0; i<3; i++)     error[i] = x_new[i] - x[i];
    

        actualPrecisionValue = errorNorm(error);
        x = x_new;
        iteration++;
    }
    //std::cout << iteration << std::endl;
    return x;

}
void printResult(std::vector<double>wynik) {
    for(int i=0; i<wynik.size(); i++) {
       std::cout << std::fixed << std::setprecision(10) << (i*h) << " " << wynik[i] << std::endl;
    }
}

int main(int argc, char * argv[]) {
    double lower[N-1], main[N], upper[N-1], b[N];
    int arg = atoi(argv[1]);
    double relaxationA[3][3]{{4.0, -1.0, 0.0}, {-1.0, 4.0, -1.0}, {0.0, -1.0, 4.0}};
    double relaxationB[3]{2,6,2};
    std::vector<double>wynik;
    
    initMatrix(lower, main, upper, b);

    if(arg == 1) {
        wynik = relaxation_Richardson_Method(relaxationA, relaxationB, 1e-10);
        for(int i=0; i<wynik.size(); i++) {
            std::cout<< (i) << std::fixed << std::setprecision(10)  << " " << wynik[i] << std::endl;
        }
    } else {
        if(arg == 2) wynik = jacobi_Method(lower, main, upper, b, 1e-10);
        if(arg == 3) wynik = gauss_Seidel_Method(lower, main, upper, b, 1e-10);
        if(arg == 4) wynik = successive_OverRelaxation_Method(lower, main, upper, b, 1e-10);
        printResult(wynik);
    }
    
}