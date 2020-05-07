#ifndef solver_master
#define solver_master

#include <iostream>
#include <complex>
using namespace std;

namespace solver{

    class RealVariable{

    public:

        double x_2;
        double x_1;
        double x_0;

        RealVariable(){
            this->x_1=1;
            this->x_2=0;
            this->x_0=0;
        }
        ~RealVariable() {}
		
        RealVariable(double x0, double x1, double x2) {
			this->x_2 = x2;
			this->x_1 = x1;
			this->x_0 = x0;
		}

        double getx_0() const{
         return this->x_0;
        }
        double getx_1() const{
         return this->x_1;
        }
        double getx_2() const{
         return this->x_2;
        }

        const RealVariable operator==(const double oth) const;
        const RealVariable operator*(const double oth) const;
        const RealVariable operator-(const double oth) const;
        const RealVariable operator+(const double oth) const;
        const RealVariable operator/(const double oth) const;
        const RealVariable operator^(const double oth) const;
        
		const RealVariable operator==(const RealVariable& re) const;
        const RealVariable operator*(const RealVariable& re) const;
        const RealVariable operator-(const RealVariable& re) const;
        const RealVariable operator+(const RealVariable& re) const;
        const RealVariable operator/(const RealVariable& re) const;

        friend const RealVariable operator==(double x,const RealVariable &im);
        friend const RealVariable operator/(double x,const RealVariable &im);
        friend const RealVariable operator+(double x,const RealVariable &im);
        friend const RealVariable operator*(double x,const RealVariable &im);
        friend const RealVariable operator-(double x,const RealVariable &im);
          
    };

    class ComplexVariable {

        std::complex<double> x_0;
        std::complex<double> x_1;
        std::complex<double> x_2;
    
    public:

       ComplexVariable(){
           this->x_0=0.0;
           this->x_1=1.0;
           this->x_2=0.0;
       }
       
       ~ComplexVariable(){}  
       
       std::complex<double> getx_0() const{
         return this->x_0;
        }
        std::complex<double> getx_1() const{
         return this->x_1;
        }
        std::complex<double> getx_2() const{
         return this->x_2;
        }
        
        const ComplexVariable operator+(std::complex<double> x) const;
        const ComplexVariable operator==(std::complex<double> x) const;
        const ComplexVariable operator-(std::complex<double> x) const;
        const ComplexVariable operator^(std::complex<double> x) const;
        const ComplexVariable operator/(std::complex<double> x) const;
        const ComplexVariable operator*(std::complex<double> x) const;

        const ComplexVariable operator*(const ComplexVariable &im) const;
        const ComplexVariable operator+(const ComplexVariable &im) const;
        const ComplexVariable operator==(const ComplexVariable &im) const;
        const ComplexVariable operator/(const ComplexVariable &im) const;
        const ComplexVariable operator-(const ComplexVariable &im) const;

        friend const ComplexVariable operator+(std::complex<double> x,const ComplexVariable &im);
        friend const ComplexVariable operator*(std::complex<double> x,const ComplexVariable &im);
        friend const ComplexVariable operator==(std::complex<double> x,const ComplexVariable &im);
        friend const ComplexVariable operator-(std::complex<double> x,const ComplexVariable &im);
        friend const ComplexVariable operator/(std::complex<double> x,const ComplexVariable &im);
    };
   
        double solve(const RealVariable& im);
        std::complex<double> solve(const ComplexVariable& im);
    }

#endif 
