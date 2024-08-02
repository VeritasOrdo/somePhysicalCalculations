#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<fstream>
#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../package/BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
#include "../package/RadiationOfElectron/BasicRadiation/BasicRadiation.h"

int main(){
    double electronMass = 511000;
    double fieldParameter1 = 20;
    double fieldParameter2 = 0.3;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 78.89*electronMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double omega = 1.55;
    std::fstream file;
    file.open("test.txt",std::ios::out);
    for(double phonEnergyRate=0.0;phonEnergyRate<0.001;phonEnergyRate+=0.00001){
        if(phonEnergyRate==0.0){
            file << phonEnergyRate << "\t" << 0 << std::endl;
            continue;
        }
        double photonEnergy = phonEnergyRate*energyPrime;
        //ElectronInCounterpropagatingLaser electronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0);
        BasicRadiationOfElectronInCounterpropagatingLaser basicRadiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,0);
        basicRadiationOfElectron.calculateDifferentialEmissionIntensity();
        //file << photonEnergy << "\t" << basicRadiationOfElectron.getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle() << std::endl;
        file << phonEnergyRate << "\t" << basicRadiationOfElectron.getDifferentialEmissionIntensity() << std::endl;
    }
    file.close();
    /*std::cout<<energyPrime<<std::endl;
    double photonEnergy = 4088*3;
    std::cout<<"photonEnengy: "<<photonEnergy<<std::endl;
    std::cout<<"momentumZPrime: "<<momentumZPrime<<std::endl;
    std::cout<<"momentumXPrime: "<<momentumXPrime<<std::endl;
    std::cout<<"fieldParameter1: "<<fieldParameter1<<std::endl;
    std::cout<<"fieldParameter2: "<<fieldParameter2<<std::endl;
    ElectronInCounterpropagatingLaser electronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0);
    std::cout<<"omega2/omega1: "<<electronInCounterpropagatingLaser.getOmega2()/electronInCounterpropagatingLaser.getOmega1()<<std::endl;
    BasicRadiationOfElectronInCounterpropagatingLaser basicRadiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,0);
    //basicRadiationOfElectron.test();
    basicRadiationOfElectron.calculateDifferentialEmissionIntensity();
    std::cout << "differentialEmissionIntensity: " << basicRadiationOfElectron.getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle() << std::endl;*/
    return 0;
}

