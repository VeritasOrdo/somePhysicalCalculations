#include "../../package/Dimension3Vector/Dimension3Vector.h"
#include<cmath>
#include<vector>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#pragma once

class NumericalSolutionOfElectronInCounterpropagatingLaser {
    protected:
        const double omega=1.55;
        const double electronMass = 511000;
    private:
        double stepLength;
        std::vector<Dimension3Vector<double>> electronMomentum;
        Dimension3Vector<double> electronMomentumPrime;
        std::vector<Dimension3Vector<double>> electronCoordinate;
        Dimension3Vector<double> electronCoordinatePrime;
        std::vector<Dimension3Vector<double>> electronVelocity;
        Dimension3Vector<double> electronVelocityPrime;
        double fieldParameter1;
        double fieldParameter2;
        double properTimeBegin;
        double properTimeEnd;
        std::vector<double> properTime;
        std::vector<double> Enengy;
        double EnengyPrime;
    public:
        NumericalSolutionOfElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTimeBegin,double properTimeEnd,double stepLength);
        std::vector<Dimension3Vector<double>> getElectronMomentum();
        Dimension3Vector<double> getElectronMomentumPrime();
        std::vector<Dimension3Vector<double>> getElectronCoordinate();
        Dimension3Vector<double> getElectronCoordinatePrime();
        std::vector<Dimension3Vector<double>> getElectronVelocity();
        Dimension3Vector<double> getElectronVelocityPrime();
        std::vector<double> getProperTime();
        std::vector<double> getEnengy();
        double getEnengyPrime();
        double getStepLength();
        double getFieldParameter1();
        double getFieldParameter2();
        double getProperTimeBegin();
        double getProperTimeEnd();
        double getOmega();
        double getElectronMass();
};