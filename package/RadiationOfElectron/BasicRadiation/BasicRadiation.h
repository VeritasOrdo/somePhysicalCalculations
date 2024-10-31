#include <cmath>
#include <iostream>
#include <complex>
#include <utility>
#include <map>

#include "../../BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
#include "../../ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#pragma once

class BasicRadiationOfElectronInCounterpropagatingLaser: public ElectronInCounterpropagatingLaser {
    private:
        double photonEnergy;
        double emissionAzimuthalAngle;
        double differentialEmissionIntensity;
        double residualEnergy;
        double energyRatio;
        double rotationDirection1;
        double rotationDirection2;
        double rotationDirectionPlus;
        double rotationDirectionMinus;
        std::vector<std::map<std::pair<int,int>,std::pair<double,double>>> emissionMapIntensityList;
        double calculateZ1X(double emissionPolarAngle);
        double calculateZ1Y(double emissionPolarAngle);
        double calculateZ2X(double emissionPolarAngle);
        double calculateZ2Y(double emissionPolarAngle);
        double trigonometricCoefficient1(double emissionPolarAngle);
        double trigonometricCoefficient2(double emissionPolarAngle);
        double trigonometricCoefficient3(double emissionPolarAngle);
        double auxiliaryAngle1(double emissionPolarAngle);
        double auxiliaryAngle2(double emissionPolarAngle);
    protected:
        
    public:
        void setDifferentialEmissionIntensity(double differentialEmissionIntensity);
        BasicRadiationOfElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime,double photonEnergy,double emissionAzimuthalAngle,double rotationDirection1,double rotationDirection2);
        virtual void calculateDifferentialEmissionIntensity();
        virtual void calculateEmissionIntensityAndPolarization();
        double getPhotonEnergy();
        double getEmissionAzimuthalAngle();
        double getDifferentialEmissionIntensity();
        double getResidualEnergy();
        double getEnergyRatio();
        std::vector<double> fourAmplitudesOfDifferentialEmissionIntensity();
        std::vector<std::map<std::pair<int,int>,std::pair<double,double>>> getEmissionMapIntensityList();
        std::vector<int> calculateLabelLimits();
        std::vector<double> calculateEmissionPolarAngle(int labelLeft, int labelRight);
        std::vector<std::complex<double>> SpectralComponent(int labelLeft, int labelRight, int label3, double emissionPolarAngle);
        void test();
        ~BasicRadiationOfElectronInCounterpropagatingLaser();
};

