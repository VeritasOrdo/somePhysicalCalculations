#include "BasicMathFunctionDefinition.h"
#include <cmath>
#include <complex>

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(double lable,double argument,double phase) {
    if(std::cyl_bessel_j(lable,argument)<1e-10) {
        return std::complex<double>(-99999999,-99999999);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::cyl_bessel_j(lable,argument);
}

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(double lable,double argument,double phase) {
    if(std::cyl_bessel_j(lable,argument)<1e-10) {
        return std::complex<double>(-99999999,-99999999);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((lable/argument)*std::cyl_bessel_j(lable,argument)*std::cos(phase),-std::sin(phase)*((1/2)*(std::cyl_bessel_j(lable-1,argument)-std::cyl_bessel_j(lable+1,argument))));
}

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(double lable,double argument,double phase) {
    if(std::cyl_bessel_j(lable,argument)<1e-10) {
        return std::complex<double>(-99999999,-99999999);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((lable/argument)*std::cyl_bessel_j(lable,argument)*std::sin(phase),std::cos(phase)*((1/2)*(std::cyl_bessel_j(lable-1,argument)-std::cyl_bessel_j(lable+1,argument))));
}
