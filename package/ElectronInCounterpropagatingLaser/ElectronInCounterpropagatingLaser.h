#include "../LorentzVector/LorentzVector.h"
#include<cmath>
#include<vector>
#pragma once

class ElectronInCounterpropagatingLaser {
    protected:
        const double omega=1.55;
        const double electronMass = 511000;
    private:
        LorentzVector *electronLorentzMomentum;
        LorentzVector *electronLorentzMomentumPrime;
        LorentzVector *electronLorentzCoordinate;
        std::vector<double> electronVelocity;
        std::vector<double> electronVelocityPrime;
        double fieldParameter1;
        double fieldParameter2;
        double omega1;
        double omega2;
        double phase1;
        double phase2;
        double properTime;
    public:
        ElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime);
        LorentzVector getElectronLorentzMomentum();
        LorentzVector getElectronLorentzMomentumPrime();
        LorentzVector getElectronLorentzCoordinate();
        std::vector<double> getElectronVelocity();
        std::vector<double> getElectronVelocityPrime();
        double getFieldParameter1();
        double getFieldParameter2();
        double getPhase1();
        double getPhase2();
        double getOmega1();
        double getOmega2();
        double getEnergy();
        double getEnergyPrime();
        double getInitialMomentumX();
        double getInitialMomentumZ();
        double getVelocityZPrime();
        double getVelocityXPrime();
        ~ElectronInCounterpropagatingLaser();
};