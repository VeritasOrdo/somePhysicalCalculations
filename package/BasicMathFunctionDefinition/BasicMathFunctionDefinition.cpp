#include "BasicMathFunctionDefinition.h"
#include <cmath>
#include <complex>
#include <iostream> 

/*std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(int lable,double argument,double phase) {
    if(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)<1e-300||std::isnan(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument))) {
        return std::complex<double>(0,0);
    }else {
        return std::exp(std::complex<double>(0,lable*phase))*BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument);
    }
}*/

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(double bessel,int label,double phase) {
    return std::exp(std::complex<double>(0,label*phase))*bessel;
}

/*std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(int lable,double argument,double phase) {
    if(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)<1e-300||std::isnan(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument))) {
        return std::complex<double>(0,0);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((lable/argument)*BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)*std::cos(phase),-std::sin(phase)*((1/2)*(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable-1,argument)-BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable+1,argument))));
}*/

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(double bessel,double besselMinus,double besselPlus,int label,double argument,double phase) {
    return std::exp(std::complex<double>(0,label*phase))*std::complex<double>((0.5)*(besselMinus+besselPlus)*std::cos(phase),-std::sin(phase)*((0.5)*(besselMinus-besselPlus)));
}

/*std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(int lable,double argument,double phase) {
    if(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)<1e-300||std::isnan(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument))) {
        return std::complex<double>(0,0);
    }
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((lable/argument)*BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable,argument)*std::sin(phase),std::cos(phase)*((1/2)*(BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable-1,argument)-BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(lable+1,argument))));
}*/

std::complex<double> BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(double bessel,double besselMinus,double besselPlus,int lable,double argument,double phase) {
    //std::cout<<std::exp(std::complex<double>(0,lable*phase))<<'\t'<<std::complex<double>((lable/argument)*bessel*std::sin(phase),std::cos(phase)*((1/2)*(besselMinus-besselPlus)))<<std::endl;
    return std::exp(std::complex<double>(0,lable*phase))*std::complex<double>((0.5)*(besselMinus+besselPlus)*std::sin(phase),std::cos(phase)*((0.5)*(besselMinus-besselPlus)));
}