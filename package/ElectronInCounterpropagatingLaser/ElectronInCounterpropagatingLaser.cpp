#include "ElectronInCounterpropagatingLaser.h"
#include<iostream>

ElectronInCounterpropagatingLaser::ElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime) {
    double momentumYPrime = 0;
    double electronMass = 511000;
    this->fieldParameter1 = fieldParameter1;
    this->fieldParameter2 = fieldParameter2;
    this->properTime = properTime;
    double reducedMass = electronMass*std::sqrt(1+this->fieldParameter1*this->fieldParameter1+this->fieldParameter2*this->fieldParameter2);
    double enengyPrime = std::sqrt(momentumXPrime*momentumXPrime+momentumYPrime*momentumYPrime+momentumZPrime*momentumZPrime+reducedMass*reducedMass);
    this->electronLorentzMomentumPrime = new LorentzVector(enengyPrime,momentumXPrime,momentumYPrime,momentumZPrime);
    this->electronVelocityPrime={(*this->electronLorentzMomentumPrime)[1]/(*this->electronLorentzMomentumPrime)[0],(*this->electronLorentzMomentumPrime)[2]/(*this->electronLorentzMomentumPrime)[0],(*this->electronLorentzMomentumPrime)[3]/(*this->electronLorentzMomentumPrime)[0]};
    this->phase1 = ((omega*(*this->electronLorentzMomentumPrime)[0]*(1-this->electronVelocityPrime[2]))/electronMass)*this->properTime;
    this->phase2 = ((omega*(*this->electronLorentzMomentumPrime)[0]*(1+this->electronVelocityPrime[2]))/electronMass)*this->properTime;
    this->electronLorentzMomentum = new LorentzVector(
        (*this->electronLorentzMomentumPrime)[0]+
        (*this->electronLorentzMomentumPrime)[1]*omega*(
            ((electronMass*fieldParameter1)/(omega*(*this->electronLorentzMomentumPrime)[0]*(1-this->electronVelocityPrime[2])))*std::cos(this->phase1)+
            ((electronMass*fieldParameter2)/(omega*(*this->electronLorentzMomentumPrime)[0]*(1+this->electronVelocityPrime[2])))*std::cos(this->phase2)
        )
        ,
        (*this->electronLorentzMomentumPrime)[1]+
        electronMass*fieldParameter1*std::cos(this->phase1)+electronMass*fieldParameter2*std::cos(this->phase2)
        ,
        (*this->electronLorentzMomentumPrime)[2]+
        electronMass*fieldParameter1*std::sin(this->phase1)+electronMass*fieldParameter2*std::sin(this->phase2)
        ,
        (*this->electronLorentzMomentumPrime)[3]+
        (*this->electronLorentzMomentumPrime)[1]*omega*(
            ((electronMass*fieldParameter1)/(omega*(*this->electronLorentzMomentumPrime)[0]*(1-this->electronVelocityPrime[2])))*std::cos(this->phase1)-
            ((electronMass*fieldParameter2)/(omega*(*this->electronLorentzMomentumPrime)[0]*(1+this->electronVelocityPrime[2])))*std::cos(this->phase2)
        )+
        ((2*electronMass*electronMass*fieldParameter1*fieldParameter2*omega)/(-2*omega*(*this->electronLorentzMomentumPrime)[0]*this->electronVelocityPrime[2]))*std::cos(this->phase1-this->phase2)
    );
    this->electronVelocity={(*this->electronLorentzMomentum)[1]/(*this->electronLorentzMomentum)[0],(*this->electronLorentzMomentum)[2]/(*this->electronLorentzMomentum)[0],(*this->electronLorentzMomentum)[3]/(*this->electronLorentzMomentum)[0]};
    this->electronLorentzCoordinate = new LorentzVector(
        ((*this->electronLorentzMomentumPrime)[0]/electronMass)*this->properTime+
        (*this->electronLorentzMomentumPrime)[1]*omega*(
            ((electronMass*fieldParameter1)/((omega*(*this->electronLorentzMomentumPrime)[0]*(1-this->electronVelocityPrime[2]))*(omega*(*this->electronLorentzMomentumPrime)[0]*(1-this->electronVelocityPrime[2]))))*std::sin(this->phase1)+
            ((electronMass*fieldParameter2)/((omega*(*this->electronLorentzMomentumPrime)[0]*(1+this->electronVelocityPrime[2]))*(omega*(*this->electronLorentzMomentumPrime)[0]*(1+this->electronVelocityPrime[2]))))*std::sin(this->phase2)
        )
        ,
        ((*this->electronLorentzMomentumPrime)[1]/electronMass)*this->properTime+
        ((electronMass*fieldParameter1)/(omega*((*this->electronLorentzMomentumPrime)[0])*(1-this->electronVelocityPrime[2])))*std::sin(this->phase1)+
        ((electronMass*fieldParameter1)/(omega*((*this->electronLorentzMomentumPrime)[0])*(1+this->electronVelocityPrime[2])))*std::sin(this->phase2)
        ,
        -
        ((electronMass*fieldParameter1)/(omega*(*this->electronLorentzMomentumPrime)[0])*(1-this->electronVelocityPrime[2]))*std::cos(this->phase1)-
        ((electronMass*fieldParameter2)/(omega*(*this->electronLorentzMomentumPrime)[0])*(1+this->electronVelocityPrime[2]))*std::cos(this->phase2)
        ,
        (((*this->electronLorentzMomentumPrime)[3])/electronMass)*this->properTime+
        ((2*electronMass*electronMass*fieldParameter1*fieldParameter2*omega)/((-2*omega*(*this->electronLorentzMomentumPrime)[0]*this->electronVelocityPrime[2])*(-2*omega*(*this->electronLorentzMomentumPrime)[0]*this->electronVelocityPrime[2])))*std::sin(this->phase1-this->phase2)+
        ((*this->electronLorentzMomentumPrime)[1])*omega*(
            ((electronMass*fieldParameter1)/((omega*(*this->electronLorentzMomentumPrime)[0]*(1-this->electronVelocityPrime[2]))*(omega*((*this->electronLorentzMomentumPrime)[0])*(1-this->electronVelocityPrime[2]))))*std::sin(this->phase1)-
            ((electronMass*fieldParameter2)/((omega*(*this->electronLorentzMomentumPrime)[0]*(1+this->electronVelocityPrime[2]))*(omega*((*this->electronLorentzMomentumPrime)[0])*(1+this->electronVelocityPrime[2]))))*std::sin(this->phase2)
        )
    );
    this->omega1 = omega*(1-this->electronVelocityPrime[2]);
    this->omega2 = omega*(1+this->electronVelocityPrime[2]);
}

LorentzVector ElectronInCounterpropagatingLaser::getElectronLorentzMomentum() {
    return *this->electronLorentzMomentum;
}

LorentzVector ElectronInCounterpropagatingLaser::getElectronLorentzMomentumPrime() {
    return *this->electronLorentzMomentumPrime;
}

LorentzVector ElectronInCounterpropagatingLaser::getElectronLorentzCoordinate() {
    return *this->electronLorentzCoordinate;
}

double ElectronInCounterpropagatingLaser::getFieldParameter1() {
    return this->fieldParameter1;
}

double ElectronInCounterpropagatingLaser::getFieldParameter2() {
    return this->fieldParameter2;
}

double ElectronInCounterpropagatingLaser::getPhase1() {
    return this->phase1;
}

double ElectronInCounterpropagatingLaser::getPhase2() {
    return this->phase2;
}

double ElectronInCounterpropagatingLaser::getOmega1() {
    return this->omega1;
}

double ElectronInCounterpropagatingLaser::getOmega2() {
    return this->omega2;
}

std::vector<double> ElectronInCounterpropagatingLaser::getElectronVelocity() {
    return this->electronVelocity;
}

std::vector<double> ElectronInCounterpropagatingLaser::getElectronVelocityPrime() {
    return this->electronVelocityPrime;
}

double ElectronInCounterpropagatingLaser::getEnergy() {
    return (*this->electronLorentzMomentum)[0];
}

double ElectronInCounterpropagatingLaser::getEnergyPrime() {
    return (*this->electronLorentzMomentumPrime)[0];
}

double ElectronInCounterpropagatingLaser::getInitialMomentumX() {
    return (*this->electronLorentzMomentumPrime)[1];
}

double ElectronInCounterpropagatingLaser::getInitialMomentumZ() {
    return (*this->electronLorentzMomentumPrime)[3];
}

double ElectronInCounterpropagatingLaser::getVelocityZPrime() {
    return this->electronVelocityPrime[2];
}

double ElectronInCounterpropagatingLaser::getVelocityXPrime() {
    return this->electronVelocityPrime[0];
}

ElectronInCounterpropagatingLaser::~ElectronInCounterpropagatingLaser() {
    delete this->electronLorentzMomentum;
    delete this->electronLorentzMomentumPrime;
    delete this->electronLorentzCoordinate;
}