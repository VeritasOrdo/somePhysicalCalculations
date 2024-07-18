#include <cmath>
#include <iostream>
#include <complex>
#include <utility>
#include "../../BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
#include "../../ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#pragma once

class BasicRadiationOfElectronInCounterpropagatingLaser: private ElectronInCounterpropagatingLaser {
    private:
        double photonEnergy;
        double emissionAzimuthalAngle;
        double differentialEmissionIntensity;
        double residualEnergy;
        double energyRatio;
        std::vector<double> calculateEmissionPolarAngle(int labelLeft, int labelRight);
        double calculateZ1X(double emissionPolarAngle);
        double calculateZ1Y(double emissionPolarAngle);
        double calculateZ2X(double emissionPolarAngle);
        double calculateZ2Y(double emissionPolarAngle);
        double trigonometricCoefficient1(double emissionPolarAngle);
        double trigonometricCoefficient2(double emissionPolarAngle);
        double trigonometricCoefficient3(double emissionPolarAngle);
        double auxiliaryAngle1(double emissionPolarAngle);
        double auxiliaryAngle2(double emissionPolarAngle);
        std::complex<double> spectralComponentT(int labelLeft, int labelRight, int label3, double emissionPolarAngle);
        std::complex<double> spectralComponentX(int labelLeft, int labelRight, int label3, double emissionPolarAngle);
        std::complex<double> spectralComponentY(int labelLeft, int labelRight, int label3, double emissionPolarAngle);
        std::complex<double> spectralComponentZ(int labelLeft, int labelRight, int label3, double emissionPolarAngle);
    public:
        BasicRadiationOfElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime,double photonEnergy,double emissionAzimuthalAngle);
        void calculateDifferentialEmissionIntensity();
        double getPhotonEnergy();
        double getEmissionAzimuthalAngle();
        double getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle();
        void test();
        ~BasicRadiationOfElectronInCounterpropagatingLaser();
};

