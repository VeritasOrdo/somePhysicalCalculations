#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<fstream>

int main() {
    double electronMass = 0.511;
    double fieldParameter1 = 50;
    double fieldParameter2 = 1;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 182*electronMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    std::fstream file;
    file.open("output.txt",std::ios::out);
    double PI = 3.14159265358979323846;
    double omega = 1;
    //std::cout << (2*PI)/(omega*(1-momentumZPrime/EnergyPrime)) << std::endl;
    for (double properTime = 0; properTime <((2*PI*electronMass)/(omega*(1-momentumZPrime/energyPrime)*energyPrime)); properTime+=((2*PI*electronMass)/(omega*(1-momentumZPrime/energyPrime)*energyPrime))/10000){
        //std::cout << properTime << std::endl;
        ElectronInCounterpropagatingLaser electron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,properTime);
        //out put the velocity of the electron
        file<<electron.getElectronVelocity()[0]<<"\t"<<electron.getElectronVelocity()[1]<<"\t"<<electron.getElectronVelocity()[2]<<" "<<std::endl;
    }
    file.close();
    return 0;
}