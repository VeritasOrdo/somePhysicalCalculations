#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<fstream>
#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../package/RadiationOfElectron/BasicRadiation/BasicRadiation.h"

int main(){
    double electronMass = 511000;
    double fieldParameter1 = 50;
    double fieldParameter2 = 0.3;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 4.0*reducedMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double photonEnergy = 0.00005*energyPrime;
    double omega = 1.55;
    std::cout<<"m*/pz: "<<reducedMass/momentumZPrime<<std::endl;
    std::cout<<"vz: "<<momentumZPrime/energyPrime<<std::endl;
    std::cout<<"omega1/omega2: "<<(1+momentumZPrime/energyPrime)/(1-momentumZPrime/energyPrime)<<std::endl;
    ElectronInCounterpropagatingLaser electronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0);
    //BasicRadiationOfElectronInCounterpropagatingLaser basicRadiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,PI/2);
    //std::cout<<"omega: "<<electronInCounterpropagatingLaser.getOmega()<<std::endl;
    std::cout<<"omega1: "<<electronInCounterpropagatingLaser.getOmega1()<<std::endl;
    std::cout<<"omega2: "<<electronInCounterpropagatingLaser.getOmega2()<<std::endl;
    std::cout<<"omega1/omega2: "<<electronInCounterpropagatingLaser.getOmega2()/electronInCounterpropagatingLaser.getOmega1()<<std::endl;
    //basicRadiationOfElectron.test();
    //basicRadiationOfElectron.calculateDifferentialEmissionIntensity();
    //std::cout << "differentialEmissionIntensity: " << basicRadiationOfElectron.getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle() << std::endl;
    return 0;
}

