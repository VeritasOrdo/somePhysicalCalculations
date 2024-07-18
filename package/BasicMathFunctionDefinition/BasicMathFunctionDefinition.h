#include<cmath>
#include<complex>
#include<utility>
#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_errno.h>
#pragma once

class BasicMathFunctionDefinition {
    public:
        static double BesselFunctionForIntegerOrder(int lable,double argument);
        static std::complex<double> relatedBesselFunctionZeroKind(int lable,double argument,double phase);
        static std::complex<double> relatedBesselFunctionFirstKind(int lable,double argument,double phase);
        static std::complex<double> relatedBesselFunctionSecondKind(int lable,double argument,double phase);
};