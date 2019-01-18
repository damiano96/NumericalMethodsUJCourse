#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <map>

int main(int argc, char const *argv[]) {
    std::map<double, std::vector<double>> mapWithValues;
    std::vector<double> values;

    double r = 2.0;
    double x = 0.5;
    int end_step = 500;
       
    while(r<=4) {
        values.clear();
        for(int j=100; j<end_step; j++) {
            x = (r*x*(1-x));
            values.push_back(x);
        }
        mapWithValues[r] = values;
        r+=0.005;
    }

    for(auto const& iterator1: mapWithValues) {
        for(auto const& iterator2: iterator1.second) {
            std::cout << std::setprecision(5) << iterator1.first << "\t" << std::setprecision(15)<< iterator2 << std::endl;
        }
    }

    return 0;
}