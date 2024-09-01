#include "../NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum/NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum.h"
#include <fstream>
#include <ctime>
#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>

/*int main(){
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
}*/

int main(){
    double omega = 1.55;
    double fieldParameter1 = 50;
    double fieldParameter2 = 1;
    double electronMass = 511000;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyInit = 250*electronMass;
    std::cout<<"energyInit: "<<energyInit<<std::endl;
    double momentumXInit = ((electronMass*fieldParameter1)+(electronMass*fieldParameter2));
    double momentumZInit = std::sqrt((energyInit*energyInit)-(electronMass*electronMass)-(momentumXInit*momentumXInit));
    std::cout<<"momentumZInit: "<<momentumZInit<<std::endl;
    double momentumZPrime = std::sqrt((energyInit*energyInit)-(reducedMass*reducedMass));
    std::cout<<"momentumZPrime : "<<momentumZPrime <<std::endl;
    double velocityXInit = momentumXInit/energyInit;
    double velocityZInit = momentumZInit/energyInit;
    double velocityZPrime = momentumZPrime/energyInit;
    double timeInterval = (M_PI*2)/(omega*(1-velocityZPrime));
    std::cout<<"velocityXInit: "<<velocityXInit<<std::endl;
    std::cout<<"velocityZInit: "<<velocityZInit<<std::endl;
    //time start
    auto start = std::chrono::high_resolution_clock::now();
    //numerical solution of molecule in counterpropagating laser*
    NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum numericalSolutionOfMoleculeInCounterpropagatingLaser(electronMass, Dimension3Vector<double>(0,0,0),Dimension3Vector<double>(momentumXInit,0,momentumZInit),0,timeInterval,timeInterval/20000.0,[omega,fieldParameter1,fieldParameter2,electronMass](double t,Dimension3Vector<double> y,Dimension3Vector<double> dydt)->Dimension3Vector<double>{  
        Dimension3Vector<double> relatedElectricField1(-fieldParameter1*std::sin(omega*t-omega*y.getZ()),fieldParameter1*std::cos(omega*t-omega*y.getZ()),0);
        Dimension3Vector<double> relatedElectricField2(-fieldParameter2*std::sin(omega*t+omega*y.getZ()),fieldParameter2*std::cos(omega*t+omega*y.getZ()),0);
        Dimension3Vector<double> relatedMagneticField1(-fieldParameter1*std::cos(omega*t-omega*y.getZ()),-fieldParameter1*std::sin(omega*t-omega*y.getZ()),0);
        Dimension3Vector<double> relatedMagneticField2(fieldParameter2*std::cos(omega*t+omega*y.getZ()),fieldParameter2*std::sin(omega*t+omega*y.getZ()),0);
        Dimension3Vector<double> relatedElectricField = relatedElectricField1+relatedElectricField2;
        Dimension3Vector<double> relatedMagneticField = relatedMagneticField1+relatedMagneticField2;
        double energy = std::sqrt((electronMass*electronMass)+(dydt*dydt));
        Dimension3Vector<double> velocity = dydt/energy;
        return ((relatedElectricField)+(velocity^relatedMagneticField))*omega*electronMass;
    });
    //time end
    auto end = std::chrono::high_resolution_clock::now();
    //output time by ms
    std::cout<<"Time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<"ms"<<std::endl;
    //std::vector<double> t = numericalSolutionOfMoleculeInCounterpropagatingLaser.getTime();
    //std::vector<Dimension3Vector<double>> position = numericalSolutionOfMoleculeInCounterpropagatingLaser.getPosition();
    //std::vector<Dimension3Vector<double>> velocity = numericalSolutionOfMoleculeInCounterpropagatingLaser.getVelocity();
    //std::vector<Dimension3Vector<double>> momentum = numericalSolutionOfMoleculeInCounterpropagatingLaser.getMomentum();
    std::vector<double> energy = numericalSolutionOfMoleculeInCounterpropagatingLaser.getEnergy();
    std::vector<Dimension3Vector<double>> velocity = numericalSolutionOfMoleculeInCounterpropagatingLaser.getVelocity();
    //file open by bin
    std::ofstream file("250mV.txt",std::ios::binary);
    //write data to file
    //time start
    auto timeFileStart = std::chrono::high_resolution_clock::now();
    //write data to file
    for(const auto& v:velocity){
        double x = v.getX();
        double y = v.getY();
        double z = v.getZ();
        file.write(reinterpret_cast<char*>(&x),sizeof(double));
        file.write(reinterpret_cast<char*>(&y),sizeof(double));
        file.write(reinterpret_cast<char*>(&z),sizeof(double));
    }
    /*for(int i=0;i<numericalSolutionOfMoleculeInCounterpropagatingLaser.getVelocity().size();i++){
        file<<numericalSolutionOfMoleculeInCounterpropagatingLaser.getVelocity()[i].getX()<<"\t"<<numericalSolutionOfMoleculeInCounterpropagatingLaser.getVelocity()[i].getY()<<"\t"<<numericalSolutionOfMoleculeInCounterpropagatingLaser.getVelocity()[i].getZ()<<std::endl;
    }*/
    auto timeFileEnd = std::chrono::high_resolution_clock::now();
    std::cout<<"TimeOfWritingFile: "<<std::chrono::duration_cast<std::chrono::milliseconds>(timeFileEnd-timeFileStart).count()<<"ms"<<std::endl;
    file.close();
    return 0;
}


