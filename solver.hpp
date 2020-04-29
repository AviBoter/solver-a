#ifndef solver_master
#define solver_master

#include <iostream>
#include <complex>
using namespace std;

namespace solver{

    class RealVariable{

    public:

        double sum = 0;
        double value;

        RealVariable& operator==(double x);

        RealVariable &operator==(RealVariable &x);

        RealVariable &operator+(double x);

        RealVariable &operator+(RealVariable &x);

        RealVariable &operator-(double x);

        RealVariable &operator-(RealVariable &x);

        RealVariable &operator*(double x);

        RealVariable &operator*(RealVariable &x);

        RealVariable &operator^(double x);

        RealVariable &operator^(RealVariable &x);

        RealVariable &operator/(double x);

        RealVariable &operator/(RealVariable &x);

    };

    class ComplexVariable {

    public:

        std::complex<double> Num = 0;
        std::complex<double> value;

        ComplexVariable &operator+(std::complex<double> x);
        ComplexVariable &operator==(std::complex<double> x);
        ComplexVariable &operator*(ComplexVariable &x);
        ComplexVariable &operator+(ComplexVariable &x);
        ComplexVariable &operator==(ComplexVariable &x);
        ComplexVariable &operator/(ComplexVariable &x);
        ComplexVariable &operator-(std::complex<double> x);
        ComplexVariable &operator^(std::complex<double> x);
        ComplexVariable &operator-(ComplexVariable &x);
        ComplexVariable &operator/(std::complex<double> x);
        ComplexVariable &operator*(std::complex<double> x);
        ComplexVariable &operator^(ComplexVariable &x); 
        
    };

    RealVariable& operator==(double x, RealVariable &y);
    RealVariable &operator/(double x, RealVariable &y);
    RealVariable &operator+(double x, RealVariable &y);
    RealVariable &operator*(double x, RealVariable &y);
    RealVariable &operator-(double x, RealVariable &y);
    
    ComplexVariable &operator+(std::complex<double> x, ComplexVariable &y);
    ComplexVariable &operator*(std::complex<double> x, ComplexVariable &y);
    ComplexVariable& operator==(std::complex<double> x, ComplexVariable &y);
    ComplexVariable &operator-(std::complex<double> x, ComplexVariable &y);
    ComplexVariable &operator/(std::complex<double> x, ComplexVariable &y);
   
     double solve(RealVariable& x);
     std::complex<double> solve(ComplexVariable& x);

};


#endif 
