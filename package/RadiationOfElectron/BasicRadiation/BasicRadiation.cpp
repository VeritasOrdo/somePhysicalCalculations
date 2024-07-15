#include "BasicRadiation.h"
#include <iostream>

BasicRadiationOfElectronInCounterpropagatingLaser::BasicRadiationOfElectronInCounterpropagatingLaser(ElectronInCounterpropagatingLaser *electronInCounterpropagatingLaser,double photonEnergy,double emissionAzimuthalAngle) {
    this->electronInCounterpropagatingLaser = electronInCounterpropagatingLaser;
    this->photonEnergy = photonEnergy;
    this->emissionAzimuthalAngle = emissionAzimuthalAngle;
    this->energyRatio = this->photonEnergy/(this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]-photonEnergy);
    this->differentialEmissionIntensity = -99999999; 
    std::cout << "Now the differential emission intensity has not been calculated." << std::endl;  
    std::cout << "Please use the function calculateDifferentialEmissionIntensity() to calculate the differential emission intensity." << std::endl;
}

std::vector<double> BasicRadiationOfElectronInCounterpropagatingLaser::calculateEmissionPolarAngle(double lableLeft,double lableRight){
    std::vector<double> emissionPolarAngle={};
    std::cout << "omega: " << omega << std::endl;
    double rho = ((this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]-photonEnergy)/(this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]*photonEnergy))*(lableLeft*(1-this->electronInCounterpropagatingLaser->getElectronVelocity()[2])*omega+lableRight*(1+this->electronInCounterpropagatingLaser->getElectronVelocity()[2])*omega);
    std::cout << "emission rho: " << rho << std::endl;
    double delta = this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*this->electronInCounterpropagatingLaser->getElectronVelocity()[2]+this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle)-(1-rho)*(1-rho);
    std::cout << "emission delta: " << delta << std::endl;
    double emissionPolarAngle1 = std::acos(
        (this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*(1-rho)+this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*std::cos(this->emissionAzimuthalAngle)*std::sqrt(delta))/
        (this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*this->electronInCounterpropagatingLaser->getElectronVelocity()[2]+this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle))
    );
    double emissionPolarAngle2 = std::acos(
        (this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*(1-rho)-this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*std::cos(this->emissionAzimuthalAngle)*std::sqrt(delta))/
        (this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*this->electronInCounterpropagatingLaser->getElectronVelocity()[2]+this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*this->electronInCounterpropagatingLaser->getElectronVelocity()[0]*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle))
    );
    std::cout << "emission polar angle 1: " << emissionPolarAngle1 << std::endl;
    std::cout << "emission polar angle 2: " << emissionPolarAngle2 << std::endl;
    if(std::sin(emissionPolarAngle1)>0){
        emissionPolarAngle.push_back(emissionPolarAngle1);
    }
    if(std::sin(emissionPolarAngle2)>0){
        emissionPolarAngle.push_back(emissionPolarAngle2);
    }
    return emissionPolarAngle;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::trigonometricCoefficient1(double lableLeft,double lableRight,double emissionPolarAngle) {
    return (electronMass*this->electronInCounterpropagatingLaser->getFieldParameter1()*energyRatio*std::sin(emissionPolarAngle))/((1-this->electronInCounterpropagatingLaser->getElectronVelocity()[2])*omega);
}

double BasicRadiationOfElectronInCounterpropagatingLaser::trigonometricCoefficient2(double lableLeft,double lableRight,double emissionPolarAngle) {
    return (electronMass*this->electronInCounterpropagatingLaser->getFieldParameter2()*energyRatio*std::sin(emissionPolarAngle))/((1+this->electronInCounterpropagatingLaser->getElectronVelocity()[2])*omega);
}

double BasicRadiationOfElectronInCounterpropagatingLaser::trigonometricCoefficient3(double lableLeft,double lableRight,double emissionPolarAngle) {
    return ((electronMass*electronMass*this->electronInCounterpropagatingLaser->getFieldParameter1()*this->electronInCounterpropagatingLaser->getFieldParameter2()*energyRatio)/(this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]*(2*this->electronInCounterpropagatingLaser->getElectronVelocity()[2]*omega)))*std::cos(this->emissionAzimuthalAngle);
}

std::complex<double> BasicRadiationOfElectronInCounterpropagatingLaser::emissionMatrixElementT(double lableLeft,double lableRight,double emissionPolarAngle) {
    std::cout << "===========================================" << std::endl;
    std::complex<double> emissionMatrixElementT = 0;
    for (int lable3=0;;lable3++){
        double lable1 = lableLeft-lable3;
        double lable2 = lableRight+lable3;
        std::cout << "lable1:"<< lable1 << std::endl;
        std::cout << "lable2:"<< lable2 << std::endl;
        std::cout << "relatedBesselFunctionZeroKind(3): " << BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable3,this->trigonometricCoefficient3(lableLeft,lableRight,emissionPolarAngle),0) << std::endl;
        std::cout << "relatedBesselFunctionZeroKind(1): " << BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle) << std::endl;
        std::cout << "relatedBesselFunctionZeroKind(2): " << BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle) << std::endl;
        std::cout << "relatedBesselFunctionFirstKind(1): " << BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle) << std::endl;
        std::cout << "relatedBesselFunctionFirstKind(2): " << BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle) << std::endl;
        if(
            BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable3,this->trigonometricCoefficient3(lableLeft,lableRight,emissionPolarAngle),0).real()<-90000000&&
            BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle).real()<-90000000&&
            BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle).real()<-90000000&&
            BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle).real()<-90000000&&
            BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle).real()<-90000000
        ){
            std::cout << "label3: " << lable3 << std::endl;
            break;
        }
        emissionMatrixElementT+=(
            BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable3,this->trigonometricCoefficient3(lableLeft,lableRight,emissionPolarAngle),0)*(
                (this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]/electronMass)*BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle)*BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle)+
                ((this->electronInCounterpropagatingLaser->getElectronLorentzMomentumPrime()[1]*omega*this->electronInCounterpropagatingLaser->getFieldParameter1())/((1-this->electronInCounterpropagatingLaser->getElectronVelocity()[2])*omega*this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]))*BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle)*BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle)+
                ((this->electronInCounterpropagatingLaser->getElectronLorentzMomentumPrime()[1]*omega*this->electronInCounterpropagatingLaser->getFieldParameter2())/((1+this->electronInCounterpropagatingLaser->getElectronVelocity()[2])*omega*this->electronInCounterpropagatingLaser->getElectronLorentzMomentum()[0]))*BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(lable1,this->trigonometricCoefficient1(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle)*BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(lable2,this->trigonometricCoefficient2(lableLeft,lableRight,emissionPolarAngle),this->emissionAzimuthalAngle)
            )
        );
    }
    return emissionMatrixElementT;
}

void BasicRadiationOfElectronInCounterpropagatingLaser::calculateDifferentialEmissionIntensity() {
    std::vector<double> emissionPolarAngle = this->calculateEmissionPolarAngle(-10,-10);
    std::cout << "emission polar angle: " << emissionPolarAngle[0] << std::endl;
    //std::cout << "===========================================" << std::endl;
    std::complex<double> emissionMatrixElementT = this->emissionMatrixElementT(0,0,emissionPolarAngle[0]);
    std::cout << "emission matrix element T: " << emissionMatrixElementT << std::endl;
    std::cout << "===========================================" << std::endl;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getPhotonEnergy() {
    return this->photonEnergy;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getEmissionAzimuthalAngle() {
    return this->emissionAzimuthalAngle;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getEmissionRelativedtoPhotonEnergyAndAzimuthalAngle() {
    return this->differentialEmissionIntensity;
}

BasicRadiationOfElectronInCounterpropagatingLaser::~BasicRadiationOfElectronInCounterpropagatingLaser() {
    delete this->electronInCounterpropagatingLaser;
}