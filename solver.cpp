
#include "solver.hpp"
#include <iostream>
#include <complex>

using namespace std;

namespace solver {

    const RealVariable RealVariable::operator==(const double oth) const {
        return RealVariable(x_0-oth,x_1,x_2);
    }

    const RealVariable RealVariable::operator*(const double oth) const {
        return RealVariable(x_0,x_1*oth,x_2);
    }
    const RealVariable RealVariable::operator-(const double x) const {
        return RealVariable(x_0-x,x_1,x_2);
    }

    const RealVariable RealVariable::operator+(const double x) const {
        return RealVariable(x_0+x,x_1,x_2);
    }

    const RealVariable RealVariable::operator^(const double x) const{
        if(x<3 && x_2==0){
            if(x_0!=0){
                return RealVariable(x_0*x_0,2*x_0*x_1,x_1*x_1);
            }
            else
            {
                return RealVariable(x_0,0,x_1);
            }   
        }
        throw std::invalid_argument("Invaild argument ^");
    }
    
     const RealVariable RealVariable::operator/(const double x) const {
        if (x==0){
            throw std::invalid_argument("double x/0");
        }
        else
        {
            return RealVariable(x_0/x,x_1/x,x_2/x);
        }
    }

    const RealVariable RealVariable::operator+(const RealVariable &x) const {
        return RealVariable(x_0+x.x_0,x_1+x.x_1,x_2+x.x_2);
    }

    const RealVariable operator==(const double x,const  RealVariable& y) {
        return y == x;
    }

    const RealVariable operator+(const double x, const RealVariable &y) {
        return y + x;
    }

    const RealVariable operator-(const double x, const RealVariable &y) {
        return x + y*-1;
    }

    const RealVariable operator*(const double x,const RealVariable &y) {
        return y * x;
    }

    const RealVariable RealVariable::operator/(const RealVariable &x) const{
        if (x.x_0!=0 && x.x_1==0 && x.x_2==0){
            return RealVariable(x_0/x.x_0,x_1/x.x_0,x_2/x.x_0);
        }
        else
        {
            throw std::invalid_argument("re x/0");
        }
    }

    ComplexVariable &ComplexVariable::operator/(ComplexVariable &x){
        return *this;

    }
    ComplexVariable& operator==(std::complex<double> x, ComplexVariable &y) {
        return y == x;
    }

    ComplexVariable& operator+(std::complex<double> x, ComplexVariable &y) {
        return y + x;
    }

    ComplexVariable& operator-(std::complex<double> x, ComplexVariable &y) {
        return x - y;
    }

    ComplexVariable& operator*(std::complex<double> x, ComplexVariable &y) {
        return y * x;
    }
    ComplexVariable& ComplexVariable::operator+(ComplexVariable &x){
        return *this;
    }
    ComplexVariable& ComplexVariable::operator==(ComplexVariable &x){
        return *this;
    }

    ComplexVariable &ComplexVariable::operator-(std::complex<double> x){
        return *this;
    }
    ComplexVariable &ComplexVariable::operator*(std::complex<double> x){
        return *this;
    }
    ComplexVariable& ComplexVariable::operator+(std::complex<double> x){
        return *this;
    }
    ComplexVariable& ComplexVariable::operator==(std::complex<double> x){
        return *this;
    }

    ComplexVariable &ComplexVariable::operator-(ComplexVariable &x){
        return *this;

    }
    ComplexVariable &ComplexVariable::operator*(ComplexVariable &x){
        return *this;

    }
     ComplexVariable &ComplexVariable::operator^(ComplexVariable &x){
        return *this;

    }
    ComplexVariable &ComplexVariable::operator^(std::complex<double> x){
        return *this;

    }
    ComplexVariable &ComplexVariable::operator/(std::complex<double> x){
        return *this;

    }
    
    double solve(RealVariable& x) {
          return 0;
    }

    std::complex<double> solve(ComplexVariable& x) {
        return x.value;
    }

};
