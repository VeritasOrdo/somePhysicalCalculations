#include "BasicMathFunctionDefinition.h"
#include <cmath>
#include <complex>
#include <iostream> 

double BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(int lable,double argument) {
    if(lable>=0||(lable%2==0)) {
        return std::cyl_bessel_j(std::abs(lable),argument);
    }else {
        return -std::cyl_bessel_j(std::abs(lable),argument);
    }
}

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(int lable,double argument,double phase) {
    if(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)<1e-300||std::isnan(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument))) {
        return std::complex<double>(0,0);
    }else {
        return std::exp(std::complex<double>(0,lable*phase))*BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument);
    }

}

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(int lable,double argument,double phase) {
    if(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)<1e-300||std::isnan(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument))) {
        return std::complex<double>(0,0);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((lable/argument)*BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)*std::cos(phase),-std::sin(phase)*((1/2)*(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable-1,argument)-BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable+1,argument))));
}

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(int lable,double argument,double phase) {
    if(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)<1e-300||std::isnan(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument))) {
        return std::complex<double>(0,0);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((lable/argument)*BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)*std::sin(phase),std::cos(phase)*((1/2)*(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable-1,argument)-BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable+1,argument))));
}
