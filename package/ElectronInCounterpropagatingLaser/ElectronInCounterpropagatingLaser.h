#include "../LorentzVector/LorentzVector.h"
#include<cmath>
#include<vector>

class ElectronInCounterpropagatingLaser {
    private:
        LorentzVector *electronLorentzMomentum;
        LorentzVector *electronLorentzMomentumPrime;
        LorentzVector *electronLorentzCoordinate;
        std::vector<double> electronVelocity;
        std::vector<double> electronVelocityPrime;
        double fieldParameter1;
        double fieldParameter2;
        double phase1;
        double phase2;
        double omega;
        double properTime;
    public:
        ElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime);
        LorentzVector getElectronLorentzMomentum();
        LorentzVector getElectronLorentzMomentumPrime();
        LorentzVector getElectronLorentzCoordinate();
        double getFieldParameter1();
        double getFieldParameter2();
        double getPhase1();
        double getPhase2();
        std::vector<double> getElectronVelocity();
        std::vector<double> getElectronVelocityPrime();
        ~ElectronInCounterpropagatingLaser();
};