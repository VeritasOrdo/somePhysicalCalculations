#include "BasicRadiation.h"
#include <iostream>
#include <algorithm>
#include <chrono>

BasicRadiationOfElectronInCounterpropagatingLaser::BasicRadiationOfElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime,double photonEnergy,double emissionAzimuthalAngle) : ElectronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,properTime) {
    this->photonEnergy = photonEnergy;
    this->emissionAzimuthalAngle = emissionAzimuthalAngle;
    this->residualEnergy = this->getEnergy()-this->photonEnergy;
    this->energyRatio = this->photonEnergy/this->residualEnergy;
    this->differentialEmissionIntensity = 0;
    std::cout << "the differential emission intensity has not been calculated yet" << std::endl;
    std::cout << "Please call the calculateDifferentialEmissionIntensity() method to calculate the differential emission intensity" << std::endl;
}

std::vector<double> BasicRadiationOfElectronInCounterpropagatingLaser::calculateEmissionPolarAngle(int labelLeft, int labelRight) {
    std::vector<double> emissionPolarAngle={};
    double rho = (1/(this->getEnergy()*energyRatio))*(labelLeft*this->getOmega1()+labelRight*this->getOmega2());
    double delta = this->getVelocityZPrime()*this->getVelocityZPrime()+this->getVelocityXPrime()*this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle)-(1-rho)*(1-rho);
    if (delta<0){
        return {};
    }
    double emissionPolarAngle1 = std::acos(
        (this->getVelocityZPrime()*(1-rho)+this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*std::sqrt(delta))/
        (this->getVelocityZPrime()*this->getVelocityZPrime()+this->getVelocityXPrime()*this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle))
    );
    double emissionPolarAngle2 = std::acos(
        (this->getVelocityZPrime()*(1-rho)-this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*std::sqrt(delta))/
        (this->getVelocityZPrime()*this->getVelocityZPrime()+this->getVelocityXPrime()*this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle))
    );
    if(std::sin(emissionPolarAngle1)>0){
        emissionPolarAngle.push_back(emissionPolarAngle1);
    }
    if(std::sin(emissionPolarAngle2)>0){
        emissionPolarAngle.push_back(emissionPolarAngle2);
    }
    return emissionPolarAngle;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::calculateZ1X(double emissionPolarAngle) {
    return ((electronMass*this->getFieldParameter1()*energyRatio)/(this->getOmega1()))*(
        std::cos(this->emissionAzimuthalAngle)*std::sin(emissionPolarAngle)+
        ((this->getVelocityXPrime()*omega)/(this->getOmega1()*this->getEnergy()))*(std::cos(emissionPolarAngle)-1)
    );
}

double BasicRadiationOfElectronInCounterpropagatingLaser::calculateZ1Y(double emissionPolarAngle) {
    return ((electronMass*this->getFieldParameter1()*energyRatio)/(this->getOmega1()))*std::sin(this->emissionAzimuthalAngle)*std::sin(emissionPolarAngle);
}

double BasicRadiationOfElectronInCounterpropagatingLaser::calculateZ2X(double emissionPolarAngle) {
    return ((electronMass*this->getFieldParameter2()*energyRatio)/(this->getOmega2()))*(
        std::cos(this->emissionAzimuthalAngle)*std::sin(emissionPolarAngle)-
        ((this->getVelocityXPrime()*omega)/(this->getOmega2()*this->getEnergy()))*(std::cos(emissionPolarAngle)+1)
    );
}

double BasicRadiationOfElectronInCounterpropagatingLaser::calculateZ2Y(double emissionPolarAngle) {
    return ((electronMass*this->getFieldParameter2()*energyRatio)/(this->getOmega2()))*std::sin(this->emissionAzimuthalAngle)*std::sin(emissionPolarAngle);
}

double BasicRadiationOfElectronInCounterpropagatingLaser::trigonometricCoefficient1(double emissionPolarAngle) {
    return std::sqrt(this->calculateZ1X(emissionPolarAngle)*this->calculateZ1X(emissionPolarAngle)+this->calculateZ1Y(emissionPolarAngle)*this->calculateZ1Y(emissionPolarAngle));
}

double BasicRadiationOfElectronInCounterpropagatingLaser::trigonometricCoefficient2(double emissionPolarAngle) {
    return std::sqrt(this->calculateZ2X(emissionPolarAngle)*this->calculateZ2X(emissionPolarAngle)+this->calculateZ2Y(emissionPolarAngle)*this->calculateZ2Y(emissionPolarAngle));
}

double BasicRadiationOfElectronInCounterpropagatingLaser::trigonometricCoefficient3(double emissionPolarAngle) {
    return ((electronMass*electronMass*this->getFieldParameter1()*this->getFieldParameter2()*energyRatio)/(this->getVelocityZPrime()*2*this->getVelocityZPrime()*omega*this->getEnergy()))*std::cos(emissionPolarAngle);
}

double BasicRadiationOfElectronInCounterpropagatingLaser::auxiliaryAngle1(double emissionPolarAngle) {
    return atan2(this->calculateZ1Y(emissionPolarAngle),this->calculateZ1X(emissionPolarAngle));
}

double BasicRadiationOfElectronInCounterpropagatingLaser::auxiliaryAngle2(double emissionPolarAngle) {
    return atan2(this->calculateZ2Y(emissionPolarAngle),this->calculateZ2X(emissionPolarAngle));
}

std::complex<double> BasicRadiationOfElectronInCounterpropagatingLaser::spectralComponentT(int labelLeft, int labelRight, int label3, double emissionPolarAngle) {
    double label1 = labelLeft-label3;
    double label2 = labelRight+label3;
    std::complex<double> zero1 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> zero2 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> zero3 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label3,this->trigonometricCoefficient3(emissionPolarAngle),0);
    std::complex<double> first1 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> first2 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> spectralComponentT = zero3*(
        (this->getEnergy()/electronMass)*zero1*zero2+
        ((this->getInitialMomentumX()*omega*this->getFieldParameter1())/(this->getOmega1()*this->getEnergy()))*first1*zero2+
        ((this->getInitialMomentumX()*omega*this->getFieldParameter2())/(this->getOmega2()*this->getEnergy()))*zero1*first2
    );
    return spectralComponentT;
}

std::complex<double> BasicRadiationOfElectronInCounterpropagatingLaser::spectralComponentX(int labelLeft, int labelRight, int label3, double emissionPolarAngle) {
    double label1 = labelLeft-label3;
    double label2 = labelRight+label3;
    std::complex<double> zero1 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> zero2 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> zero3 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label3,this->trigonometricCoefficient3(emissionPolarAngle),0);
    std::complex<double> first1 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> first2 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> spectralComponentX = zero3*(
        (this->getInitialMomentumX()/electronMass)*zero1*zero2+
        this->getFieldParameter1()*first1*zero2-
        this->getFieldParameter2()*zero1*first2
    );
    return spectralComponentX;
} 

std::complex<double> BasicRadiationOfElectronInCounterpropagatingLaser::spectralComponentY(int labelLeft, int labelRight, int label3, double emissionPolarAngle) {
    double label1 = labelLeft-label3;
    double label2 = labelRight+label3;
    std::complex<double> zero1 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> zero2 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> zero3 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label3,this->trigonometricCoefficient3(emissionPolarAngle),0);
    std::complex<double> second1 = BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> second2 = BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> spectralComponentY = zero3*(
        this->getFieldParameter1()*second1*zero2+
        this->getFieldParameter2()*zero1*second2
    );
    return spectralComponentY;
}

std::complex<double> BasicRadiationOfElectronInCounterpropagatingLaser::spectralComponentZ(int labelLeft, int labelRight, int label3, double emissionPolarAngle) {
    double label1 = labelLeft-label3;
    double label2 = labelRight+label3;
    std::complex<double> zero1 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> zero2 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> zero3 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(label3,this->trigonometricCoefficient3(emissionPolarAngle),0);
    std::complex<double> first1 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label1,this->trigonometricCoefficient1(emissionPolarAngle),this->auxiliaryAngle1(emissionPolarAngle));
    std::complex<double> first2 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label2,this->trigonometricCoefficient2(emissionPolarAngle),this->auxiliaryAngle2(emissionPolarAngle));
    std::complex<double> first3 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(label3,this->trigonometricCoefficient3(emissionPolarAngle),0);
    std::complex<double> spectralComponentZ = (
        zero1*zero2*(
            (this->getInitialMomentumZ()/electronMass)*zero3-
            ((electronMass*this->getFieldParameter1()*this->getFieldParameter2())/(this->getVelocityZPrime()*this->getEnergy()))*first3
        )+
        zero3*((this->getInitialMomentumX()*omega)/electronMass)*(
            (this->getFieldParameter1()/this->getOmega1())*first1*zero2-
            (this->getFieldParameter2()/this->getOmega2())*zero1*first2
        )
    );
    return spectralComponentZ;
}

void BasicRadiationOfElectronInCounterpropagatingLaser::calculateDifferentialEmissionIntensity() {
    double label1Limit = 0;
    double label2Limit = 0;
    double label3Limit = 0;
    for(int label1 = 0;;label1++){
        std::vector<double> trigonometricCoefficient1 = {};
        for(double polarAngle=0;polarAngle<M_PI;polarAngle+=0.01){
            trigonometricCoefficient1.push_back(this->trigonometricCoefficient1(polarAngle));
        }
        double maxTrigonometricCoefficient1 = *std::max_element(trigonometricCoefficient1.begin(),trigonometricCoefficient1.end());
        if(std::cyl_bessel_i(label1,maxTrigonometricCoefficient1)<1e-100||std::isnan(std::cyl_bessel_i(label1,maxTrigonometricCoefficient1))){
            label1Limit = label1;
            break;
        }
    }
    for(int label2 = 0;;label2++){
        std::vector<double> trigonometricCoefficient2 = {};
        for(double polarAngle=0;polarAngle<M_PI;polarAngle+=0.01){
            trigonometricCoefficient2.push_back(this->trigonometricCoefficient2(polarAngle));
        }
        double maxTrigonometricCoefficient2 = *std::max_element(trigonometricCoefficient2.begin(),trigonometricCoefficient2.end());
        if(std::cyl_bessel_i(label2,maxTrigonometricCoefficient2)<1e-100||std::isnan(std::cyl_bessel_i(label2,maxTrigonometricCoefficient2))){
            label2Limit = label2;
            break;
        }
    }
    for(int label3=0;;label3++){
        double maxTrigonometricCoefficient3 = this->trigonometricCoefficient3(M_PI/2);
        if(std::cyl_bessel_i(label3,maxTrigonometricCoefficient3)<1e-100){
            label3Limit = label3;
            break;
        }
    }
    double labelLeftLimit = label1Limit+label3Limit;
    double labelRightLimit = label2Limit+label3Limit;
    std::cout<<"labelLeftLimit: "<<labelLeftLimit<<std::endl;
    std::cout<<"labelRightLimit: "<<labelRightLimit<<std::endl;
    std::cout<<"label3Limit: "<<label3Limit<<std::endl;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    double differentialEmissionIntensityRe = 0;
    for(int labelLeft = -labelLeftLimit;labelLeft<=labelLeftLimit;labelLeft++){
        for(int labelRight = -labelRightLimit;labelRight<=labelRightLimit;labelRight++){
            std::vector<double> emissionPolarAngle = this->calculateEmissionPolarAngle(labelLeft,labelRight);
            for(int emissionPolarAngleIndex=0;emissionPolarAngleIndex<emissionPolarAngle.size();emissionPolarAngleIndex++){
                std::complex<double> spectralComponentT =0;
                std::complex<double> spectralComponentX =0;
                std::complex<double> spectralComponentY =0;
                std::complex<double> spectralComponentZ =0;
                for(int label3 = -label3Limit;label3<=label3Limit;label3++){
                    //std::cout<<labelLeft<<","<<labelRight<<","<<label3<<std::endl;
                    spectralComponentT += this->spectralComponentT(labelLeft,labelRight,label3,emissionPolarAngle[emissionPolarAngleIndex]);
                    spectralComponentX += this->spectralComponentX(labelLeft,labelRight,label3,emissionPolarAngle[emissionPolarAngleIndex]);
                    spectralComponentY += this->spectralComponentY(labelLeft,labelRight,label3,emissionPolarAngle[emissionPolarAngleIndex]);
                    spectralComponentZ += this->spectralComponentZ(labelLeft,labelRight,label3,emissionPolarAngle[emissionPolarAngleIndex]);
                }
                double spectralComponent = std::abs(spectralComponentT)*std::abs(spectralComponentT)-std::abs(spectralComponentX)*std::abs(spectralComponentX)-std::abs(spectralComponentY)*std::abs(spectralComponentY)-std::abs(spectralComponentZ)*std::abs(spectralComponentZ);
                differentialEmissionIntensityRe += spectralComponent*(1/(this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*(1/std::tan(emissionPolarAngle[emissionPolarAngleIndex]))-this->getVelocityZPrime()));
                std::cout<<"spectral component"<<spectralComponent<<std::endl;
                if(spectralComponent>1){
                    std::cout<<"spectral component"<<spectralComponent<<std::endl;
                    std::cout<<"labelLeft: "<<labelLeft<<" labelRight: "<<labelRight<<std::endl;
                }
            }
        }
    }
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    //print the time with min and sec
    std::cout<<"differentialEmissionIntensityRe: "<<differentialEmissionIntensityRe<<std::endl;
    std::cout<<"Time: "<<(time1-time0)/60000000000<<" min"<<std::endl;
    this->differentialEmissionIntensity = (((double(1)/double(137))*photonEnergy*electronMass*electronMass)/(2*M_PI*this->getEnergy()*this->getEnergy()))*differentialEmissionIntensityRe;
    //this->differentialEmissionIntensity = (((double(1)/double(137))*photonEnergy*electronMass*electronMass)/(2*M_PI*this->getEnergy()*this->getEnergy()));
    std::cout<<this->differentialEmissionIntensity<<std::endl;
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
    std::cout << "The object has been deleted" << std::endl;
}
