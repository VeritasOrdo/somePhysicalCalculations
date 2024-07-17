#include<cmath>
#include<complex>
#include<utility>
#pragma once

class BasicMathFunctionDefinition {
    public:
        static double BesselFunctionForIntegerOrder(int lable,double argument);
        static std::complex<double> relatedBesselFunctionZeroKind(int lable,double argument,double phase);
        static std::complex<double> relatedBesselFunctionFirstKind(int lable,double argument,double phase);
        static std::complex<double> relatedBesselFunctionSecondKind(int lable,double argument,double phase);
};