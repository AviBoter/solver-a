
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
    const RealVariable RealVariable::operator-(const double oth) const {
        return RealVariable(x_0-oth,x_1,x_2);
    }

    const RealVariable RealVariable::operator+(const double oth) const {
        return RealVariable(x_0+oth,x_1,x_2);
    }

    const RealVariable RealVariable::operator^(const double oth) const{
        if(oth<3 && x_2==0){
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
    
     const RealVariable RealVariable::operator/(const double oth) const {
        if (oth==0){
            throw std::invalid_argument("double x/0");
        }
        else
        {
            return RealVariable(x_0/oth,x_1/oth,x_2/oth);
        }
    }

    const RealVariable RealVariable::operator+(const RealVariable &oth) const {
        return RealVariable(x_0+oth.x_0,x_1+oth.x_1,x_2+oth.x_2);
    }

    const RealVariable operator==(const double oth,const  RealVariable& re) {
        return re == oth;
    }

    const RealVariable operator+(const double oth, const RealVariable &re) {
        return re + oth;
    }

    const RealVariable operator-(const double oth, const RealVariable &re) {
        return oth + re*-1;
    }

    const RealVariable operator*(const double oth,const RealVariable &re) {
        return re * oth;
    }

    const RealVariable RealVariable::operator/(const RealVariable &re) const{
        if (re.x_0!=0 && re.x_1==0 && re.x_2==0){
            return RealVariable(x_0/re.x_0,x_1/re.x_0,x_2/re.x_0);
        }
        else
        {
            throw std::invalid_argument("re x/0");
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     const ComplexVariable ComplexVariable::operator+(const std::complex<double> oth)const {
         ComplexVariable c;
         c.x_0=c.x_0+oth;
        return c;
    }
    const ComplexVariable ComplexVariable::operator==(const std::complex<double> oth)const {
         ComplexVariable c;
         c.x_0=c.x_0-oth;
        return c;
    }
    const ComplexVariable ComplexVariable::operator-(const std::complex<double> oth)const {
        ComplexVariable c;
        c.x_0=c.x_0-oth;
        return c;
    }

    const ComplexVariable ComplexVariable::operator*(const std::complex<double> oth)const {
         ComplexVariable c;
        c.x_0=c.x_0*oth;
        c.x_1=c.x_1*oth;
        c.x_2=c.x_2*oth;
        return c;
    }

    const ComplexVariable ComplexVariable::operator^(const std::complex<double> re)const {
        if (re==-2.0){
            throw std::invalid_argument("x^-2 isn't valid");
        }
        if(re.real()<=2.0 && x_2==0.0){
            ComplexVariable c;
            if(x_0==0.0){
              c.x_0=x_0;
              c.x_1=0;
              c.x_2=x_1;
                return c;
            }
            else
            {
              c.x_0=x_0*x_0;
              c.x_1=2.0*x_0*x_1;
              c.x_2=x_1*x_1;
                return c;
            }
        }
        throw std::invalid_argument("^ isnt valid");
    }

    const ComplexVariable ComplexVariable::operator/(const std::complex<double> re)const {
         ComplexVariable c;
        if (re==0.0){
             throw std::invalid_argument("ERRERP x/0");
            
        }
        else
        {
            c.x_0=x_0/re;
            c.x_0=x_1/re;
            c.x_0=x_2/re;
              return c;
        }

    }
  

    
    double solve(RealVariable& x) {
          return 0;
    }

    std::complex<double> solve(ComplexVariable& x) {
        return x.value;
    }

};
