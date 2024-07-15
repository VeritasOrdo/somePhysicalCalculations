#include <cmath>
#include <iostream>
#include <complex>
#include "../../ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../../BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"

class BasicRadiationOfElectronInCounterpropagatingLaser {
    private:
        const double electronMass = 0.511;
        const double omega = 1.55;
        ElectronInCounterpropagatingLaser *electronInCounterpropagatingLaser;
        double photonEnergy;
        double emissionAzimuthalAngle;
        double differentialEmissionIntensity;
        double energyRatio;
        std::vector<double> calculateEmissionPolarAngle(double lableLeft,double lableRight);
        double trigonometricCoefficient1(double lableLeft,double lableRight, double emissionPolarAngle);
        double trigonometricCoefficient2(double lableLeft,double lableRight, double emissionPolarAngle);
        double trigonometricCoefficient3(double lableLeft,double lableRight, double emissionPolarAngle);
        std::complex<double> emissionMatrixElementT(double lableLeft,double lableRight, double emissionPolarAngle);
        //std::complex<double> emissionMatrixElementX(double lableLeft,double lableRight, double emissionPolarAngle);
        //std::complex<double> emissionMatrixElementY(double lableLeft,double lableRight, double emissionPolarAngle);
        //std::complex<double> emissionMatrixElementZ(double lableLeft,double lableRight, double emissionPolarAngle);
    public:
        BasicRadiationOfElectronInCounterpropagatingLaser(ElectronInCounterpropagatingLaser *electronInCounterpropagatingLaser,double photonEnergy,double emissionAzimuthalAngle);
        void calculateDifferentialEmissionIntensity();
        double getPhotonEnergy();
        double getEmissionAzimuthalAngle();
        double getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle();
        ~BasicRadiationOfElectronInCounterpropagatingLaser();
};