#include<cmath>
#include<complex>
#include<utility>
#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_errno.h>
#pragma once

class BasicMathFunctionDefinition {
    public:
        static std::complex<double> relatedBesselFunctionZeroKind(double bessel, int lable, double phase);
        static std::complex<double> relatedBesselFunctionFirstKind(double bessel, double besselMinus, double besselPlus, int lable, double argument, double phase);
        static std::complex<double> relatedBesselFunctionSecondKind(double bessel, double besselMinus, double besselPlus, int lable, double argument, double phase);
};