#include "../package/RadiationOfElectron/BasicRadiation/BasicRadiation.h"
#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<fstream>

/*int main() {
    double electronMass = 0.511;
    double fieldParameter1 = 20;
    double fieldParameter2 = 0.3;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 4*reducedMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double photonEnergy = 100*energyPrime;
    double omega = 1.55;
    ElectronInCounterpropagatingLaser *electronInCounterpropagatingLaser = new ElectronInCounterpropagatingLaser(momentumZPrime, momentumXPrime, fieldParameter1, fieldParameter2, 0);
    BasicRadiationOfElectronInCounterpropagatingLaser *basicRadiationOfElectronInCounterpropagatingLaser = new BasicRadiationOfElectronInCounterpropagatingLaser(electronInCounterpropagatingLaser, photonEnergy, PI/2);
    basicRadiationOfElectronInCounterpropagatingLaser->calculateDifferentialEmissionIntensity();
    std::cout << "===========================================" << std::endl;
    //std::fstream file;
    //file.open("output.txt",std::ios::out);
    //file.close();
    return 0;
}*/

int main(){
    std::cout << "Bessel(0,0)" <<std::sph_bessel(-100,10) << std::endl;
}
