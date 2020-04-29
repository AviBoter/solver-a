
#include "solver.hpp"
#include <iostream>
#include <complex>

using namespace std;

namespace solver {

     
    RealVariable& operator+(double x, RealVariable &y) {
        return y + x;
    }

    RealVariable& operator==(double x, RealVariable& y) {
        return y == x;
    }
    RealVariable& operator-(double x, RealVariable &y) {
        return x + y;
    }

    RealVariable& operator*(double x, RealVariable &y) {
        return y * x;
    }

    RealVariable& RealVariable::operator+(double x) {
        return *this;
    }

    RealVariable& RealVariable::operator==(double x) {
        return *this;
    }
    RealVariable& RealVariable::operator==(RealVariable& x) {
        return *this;
    }

    RealVariable& RealVariable::operator+(RealVariable &x) {
        return *this;
    }

    RealVariable& RealVariable::operator-(RealVariable &x) {
        return *this;
    }

    RealVariable& RealVariable::operator-(double x) {
        return *this;
    }

    RealVariable& RealVariable::operator*(double x) {
        return *this ;
    }

    RealVariable& RealVariable::operator*(RealVariable &x) {
        return *this ;
    }

    RealVariable& RealVariable::operator^(double x) {
        return *this;
    }

    RealVariable& RealVariable::operator/(RealVariable &x) {
        return *this;
    }

    RealVariable& RealVariable::operator/(double x) {
        return *this;
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
