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
    double fieldParameter1 = 20;
    double fieldParameter2 = 0.3;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 4*reducedMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double photonEnergy = 0.0001*energyPrime;
    double omega = 1.55;
    //ElectronInCounterpropagatingLaser electronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0);
    BasicRadiationOfElectronInCounterpropagatingLaser basicRadiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,PI/2);
    //basicRadiationOfElectron.test();
    basicRadiationOfElectron.calculateDifferentialEmissionIntensity();
    std::cout << "differentialEmissionIntensity: " << basicRadiationOfElectron.getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle() << std::endl;
    return 0;
}

