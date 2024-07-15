#include<cmath>
#include<complex>

class BasicMathFunctionDefinition {
    public:
        static std::complex<double> relatedBesselFunctionZeroKind(double lable,double argument,double phase);
        static std::complex<double> relatedBesselFunctionFirstKind(double lable,double argument,double phase);
        static std::complex<double> relatedBesselFunctionSecondKind(double lable,double argument,double phase);
};