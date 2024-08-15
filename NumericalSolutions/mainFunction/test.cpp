#include "../NumericalSolutionOfElectronInCounterpropagatingLaser/NumericalSolutionOfElectronInCounterpropagatingLaser.h"
#include <fstream>
#include <ctime>
#include <iostream>
#include <cmath>
#include <vector>

int main(){
    double electronMass = 511000;
    double fieldParameter1 = 50;
    double fieldParameter2 = 1;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 80*reducedMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double omega = 1.55;
    double stepLength = 1e-6;
    double properTimeBegin = 0;
    double properTimeEnd = M_PI/omega;
    std::clock_t start = std::clock();
    NumericalSolutionOfElectronInCounterpropagatingLaser numericalSolutionOfElectronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,properTimeBegin,properTimeEnd,stepLength);
    std::clock_t end = std::clock();
    std::cout<<"Time: "<<(end-start)/(double)CLOCKS_PER_SEC<<"s"<<std::endl;
    std::vector<Dimension3Vector<double>> electronMomentum = numericalSolutionOfElectronInCounterpropagatingLaser.getElectronMomentum();
    std::vector<Dimension3Vector<double>> electronCoordinate = numericalSolutionOfElectronInCounterpropagatingLaser.getElectronCoordinate();
    std::vector<Dimension3Vector<double>> electronVelocity = numericalSolutionOfElectronInCounterpropagatingLaser.getElectronVelocity();
    std::vector<double> properTime = numericalSolutionOfElectronInCounterpropagatingLaser.getProperTime();
    std::vector<double> Enengy = numericalSolutionOfElectronInCounterpropagatingLaser.getEnengy();
    std::fstream file;
    file.open("electronVelocity.txt",std::ios::out);
    for(int i=0;i<electronVelocity.size();i++){
        file<<properTime[i]<<"\t"<<electronVelocity[i].getX()<<"\t"<<electronVelocity[i].getY()<<"\t"<<electronVelocity[i].getZ()<<std::endl;
    }
    file.close();
    file.open("electronCoordinate.txt",std::ios::out);
    for(int i=0;i<electronCoordinate.size();i++){
        file<<properTime[i]<<"\t"<<electronCoordinate[i].getX()<<"\t"<<electronCoordinate[i].getY()<<"\t"<<electronCoordinate[i].getZ()<<std::endl;
    }
    file.close();
    file.open("electronMomentum.txt",std::ios::out);
    for(int i=0;i<electronMomentum.size();i++){
        file<<properTime[i]<<"\t"<<electronMomentum[i].getX()<<"\t"<<electronMomentum[i].getY()<<"\t"<<electronMomentum[i].getZ()<<std::endl;
    }
    file.close();
    file.open("Enengy.txt",std::ios::out);
    for(int i=0;i<Enengy.size();i++){
        file<<properTime[i]<<"\t"<<Enengy[i]<<std::endl;
    }
    file.close();
    return 0;
}


