#include "BasicMathFunctionDefinition.h"
#include <cmath>
#include <complex>
#include <iostream> 

void my_gsl_error_handler(const char * reason, const char * file, int line, int gsl_errno) {
    //fprintf(stderr, "A GSL error occurred: %s\n", reason);
    //fprintf(stderr, "Error code: %d, in file: %s, line: %d\n", gsl_errno, file, line);
    // 根据错误类型决定是否终止程序
    if (gsl_errno != GSL_EUNDRFLW) { // 如果不是下溢错误，终止程序
        abort();
    }
    // 对于下溢错误，可以选择不终止程序
}

double BasicMathFunctionDefinition::BesselFunctionForIntegerOrder(int label,double argument) {
    gsl_set_error_handler(&my_gsl_error_handler);
    if(label>=0||(label%2==0)) {
        return gsl_sf_bessel_Jn(std::abs(label),argument);
    }else {
        return -gsl_sf_bessel_Jn(std::abs(label),argument);
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
