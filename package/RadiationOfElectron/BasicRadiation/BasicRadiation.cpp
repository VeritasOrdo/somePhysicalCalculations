#include "BasicRadiation.h"
#include "../../ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../../BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
#include <iostream>
#include <fstream>  
#include <algorithm>
#include <chrono>
#include <omp.h>


void my_gsl_error_handler(const char * reason, const char * file, int line, int gsl_errno) {
    //fprintf(stderr, "A GSL error occurred: %s\n", reason);
    //fprintf(stderr, "Error code: %d, in file: %s, line: %d\n", gsl_errno, file, line);
    // 根据错误类型决定是否终止程序
    if (gsl_errno != GSL_EUNDRFLW) { // 如果不是下溢错误，终止程序
        abort();
    }
    // 对于下溢错误，可以选择不终止程序
}


double myBesselFunction(int label,double argument) {
    gsl_set_error_handler(&my_gsl_error_handler);
    if(std::abs(argument)<1000){
        //std::cout<<"here2"<<std::endl;
        return gsl_sf_bessel_Jn(label,argument);
    }
    if(label>=0||((label%2)==0)) {
        //std::cout<<"here3"<<std::endl;
        return gsl_sf_bessel_Jnu(std::abs(label),argument);
    }else {
        //std::cout<<"here4"<<std::endl;
        return -gsl_sf_bessel_Jnu(std::abs(label),argument);
    }
}

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
    //std::cout<<"Energy: "<<this->getEnergy()<<std::endl;
    //std::cout<<"EnergyPrime: "<<this->getEnergyPrime()<<std::endl;
    //std::cout<<"energyRatio: "<<energyRatio<<"energy: "<<this->getEnergy()<<"omega1: "<<this->getOmega1()<<"omega2: "<<this->getOmega2()<<std::endl;
    //std::cout<<"labelLeft: "<<labelLeft<<" labelRight: "<<labelRight<<std::endl;
    //std::cout<<"rho: "<<rho<<std::endl;
    double delta = this->getVelocityZPrime()*this->getVelocityZPrime()+this->getVelocityXPrime()*this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*std::cos(this->emissionAzimuthalAngle)-(1-rho)*(1-rho);
    if (delta<0||(1-rho)/this->getVelocityZPrime()>1){
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
    if(std::abs(this->getVelocityXPrime())<1e-300){
        //std::cout<<"here!"<<std::endl;
        return {emissionPolarAngle1};
    }
    if(std::sin(emissionPolarAngle1)>0){
        emissionPolarAngle.push_back(emissionPolarAngle1);
    }
    if(std::sin(emissionPolarAngle2)>0){
        emissionPolarAngle.push_back(emissionPolarAngle2);
    }
    return emissionPolarAngle;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::calculateZ1X(double emissionPolarAngle) {
    //std::cout<<"========================"<<std::endl;
    //std::cout<<"electronMass: "<<electronMass<<" fieldParameter1: "<<this->getFieldParameter1()<<" energyRatio: "<<energyRatio<<" omega: "<<omega<<" getOmega1: "<<this->getOmega1()<<" getEnergy: "<<this->getEnergy()<<" getVelocityXPrime: "<<this->getVelocityXPrime()<<" emissionAzimuthalAngle: "<<this->emissionAzimuthalAngle<<" emissionPolarAngle: "<<emissionPolarAngle<<std::endl;
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
    //std::cout<<"trigonometricCoefficient3: "<<((electronMass*electronMass*this->getFieldParameter1()*this->getFieldParameter2()*energyRatio)/(this->getVelocityZPrime()*(this->getOmega2()-this->getOmega1())*this->getEnergy()))<<std::endl;
    return ((electronMass*electronMass*this->getFieldParameter1()*this->getFieldParameter2()*energyRatio)/(this->getVelocityZPrime()*(this->getOmega2()-this->getOmega1())*this->getEnergy()))*std::cos(emissionPolarAngle); 
}

double BasicRadiationOfElectronInCounterpropagatingLaser::auxiliaryAngle1(double emissionPolarAngle) {
    return atan(this->calculateZ1Y(emissionPolarAngle)/this->calculateZ1X(emissionPolarAngle));
}

double BasicRadiationOfElectronInCounterpropagatingLaser::auxiliaryAngle2(double emissionPolarAngle) {
    return atan(this->calculateZ2Y(emissionPolarAngle)/this->calculateZ2X(emissionPolarAngle));
}

std::vector<std::complex<double>> BasicRadiationOfElectronInCounterpropagatingLaser::SpectralComponent(int labelLeft, int labelRight, int label3, double emissionPolarAngle) {
    double label1 = labelLeft-label3;
    double label2 = labelRight+label3;
    double trigonometricCoefficient1 = this->trigonometricCoefficient1(emissionPolarAngle);
    double trigonometricCoefficient2 = this->trigonometricCoefficient2(emissionPolarAngle);
    double trigonometricCoefficient3 = this->trigonometricCoefficient3(emissionPolarAngle);
    //std::cout<<"label1: "<<label1<<" label2: "<<label2<<" label3: "<<label3<<" z1: "<<trigonometricCoefficient1<<" z2: "<<trigonometricCoefficient2<<" z3: "<<trigonometricCoefficient3<<std::endl;
    double auxiliaryAngle1 = this->auxiliaryAngle1(emissionPolarAngle);
    double auxiliaryAngle2 = this->auxiliaryAngle2(emissionPolarAngle);
    double bessel1 = myBesselFunction(label1,trigonometricCoefficient1);
    //std::cout<<"bessel1: "<<bessel1<<std::endl;
    double bessel2 = myBesselFunction(label2,trigonometricCoefficient2);
    //std::cout<<"bessel2: "<<bessel2<<std::endl;
    double bessel3 = myBesselFunction(label3,trigonometricCoefficient3);
    //std::cout<<"bessel3: "<<bessel3<<std::endl;
    double bessel1Plus = myBesselFunction(label1+1,trigonometricCoefficient1);
    //std::cout<<"bessel1Plus: "<<bessel1Plus<<std::endl;
    double bessel2Plus = myBesselFunction(label2+1,trigonometricCoefficient2);
    //std::cout<<"bessel2Plus: "<<bessel2Plus<<std::endl;
    double bessel3Plus = myBesselFunction(label3+1,trigonometricCoefficient3);
    //std::cout<<"bessel3Plus: "<<bessel3Plus<<std::endl;
    double bessel1Minus = myBesselFunction(label1-1,trigonometricCoefficient1);
    //std::cout<<"bessel1Minus: "<<bessel1Minus<<std::endl;
    double bessel2Minus = myBesselFunction(label2-1,trigonometricCoefficient2);
    //std::cout<<"bessel2Minus: "<<bessel2Minus<<std::endl;
    double bessel3Minus = myBesselFunction(label3-1,trigonometricCoefficient3);
    //std::cout<<"bessel3Minus: "<<bessel3Minus<<std::endl;
    std::complex<double> zero1 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(bessel1,label1,auxiliaryAngle1);
    std::complex<double> zero2 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(bessel2,label2,auxiliaryAngle2);
    std::complex<double> zero3 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(bessel3,label3,0);
    std::complex<double> first1 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(bessel1,bessel1Minus,bessel1Plus,label1,trigonometricCoefficient1,auxiliaryAngle1);
    std::complex<double> first2 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(bessel2,bessel2Minus,bessel2Plus,label2,trigonometricCoefficient2,auxiliaryAngle2);
    std::complex<double> first3 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(bessel3,bessel3Minus,bessel3Plus,label3,trigonometricCoefficient3,0);
    std::complex<double> second1 = BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(bessel1,bessel1Minus,bessel1Plus,label1,trigonometricCoefficient1,auxiliaryAngle1);
    std::complex<double> second2 = BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(bessel2,bessel2Minus,bessel2Plus,label2,trigonometricCoefficient2,auxiliaryAngle2);
    //std::cout<<"zero1: "<<zero1<<std::endl;
    //std::cout<<"zero2: "<<zero2<<std::endl;
    //std::cout<<"zero3: "<<zero3<<std::endl;
    //std::cout<<"first1: "<<first1<<std::endl;
    //std::cout<<"first2: "<<first2<<std::endl;
    //std::cout<<"first3: "<<first3<<std::endl;
    //std::cout<<"second1: "<<second1<<std::endl;
    //std::cout<<"second2: "<<second2<<std::endl;
    //std::cout<<second1<<std::endl;
    //std::cout<<second2<<std::endl;
    std::complex<double> SpectralComponentT = (
        zero3*(
            (this->getEnergy()/electronMass)*zero1*zero2+
            ((this->getInitialMomentumX()*omega*this->getFieldParameter1())/(this->getOmega1()*this->getEnergy()))*first1*zero2+
            ((this->getInitialMomentumX()*omega*this->getFieldParameter2())/(this->getOmega2()*this->getEnergy()))*zero1*first2
        )
    );
    std::complex<double> SpectralComponentX = (
        zero3*(
            (this->getInitialMomentumX()/electronMass)*zero1*zero2+
            this->getFieldParameter1()*first1*zero2+
            this->getFieldParameter2()*zero1*first2
        )
    );
    std::complex<double> SpectralComponentY = (zero3*(
            this->getFieldParameter1()*second1*zero2+
            this->getFieldParameter2()*zero1*second2
        )
    );
    std::complex<double> SpectralComponentZ = (
        zero1*zero2*(
            (this->getInitialMomentumZ()/electronMass)*zero3-
            ((electronMass*this->getFieldParameter1()*this->getFieldParameter2())/(this->getVelocityZPrime()*this->getEnergy()))*first3
        )+
        zero3*((this->getInitialMomentumX()*omega)/electronMass)*(
            (this->getFieldParameter1()/this->getOmega1())*first1*zero2-
            (this->getFieldParameter2()/this->getOmega2())*zero1*first2
        )
    );
    std::vector<std::complex<double>> spectralComponent = {SpectralComponentT,SpectralComponentX,SpectralComponentY,SpectralComponentZ};
    return spectralComponent;
}

void BasicRadiationOfElectronInCounterpropagatingLaser::calculateDifferentialEmissionIntensity() {
    double label1Limit = 0;
    double label2Limit = 0;
    double label3Limit = 0;
    double trigonometricCoefficient1Max = this->trigonometricCoefficient1(M_PI/2);
    double trigonometricCoefficient2Max = this->trigonometricCoefficient2(M_PI/2);
    double trigonometricCoefficient3Max = this->trigonometricCoefficient3(0);
    std::cout<<"z1Max: "<<trigonometricCoefficient1Max<<std::endl;
    std::cout<<"z2Max: "<<trigonometricCoefficient2Max<<std::endl;
    std::cout<<"z3Max: "<<trigonometricCoefficient3Max<<std::endl;
    long double timej0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    for(int i = 0;;i++){
        if(std::abs(myBesselFunction(i,trigonometricCoefficient1Max))<1e-320){
            label1Limit = i;
            break;
        }
    }
    for(int i = 0;;i++){
        if(std::abs(myBesselFunction(i,trigonometricCoefficient2Max))<1e-320){
            label2Limit = i;
            break;
        }
    }
    for(int i = 0;;i++){
        if(std::abs(myBesselFunction(i,trigonometricCoefficient3Max))<1e-320){
            label3Limit = i;
            break;
        }
    }
    long double timej1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"Timej: "<<(timej1-timej0)/60000000000<<" min"<<std::endl;
    int labelLeftLimit = label1Limit+label3Limit;
    int labelRightLimit = label2Limit+label3Limit;
    //labelLeftLimit = 20000;
    //labelRightLimit = 20;
    //label3Limit = 10;
    std::cout<<"labelLeftLimit: "<<labelLeftLimit<<std::endl;
    std::cout<<"labelRightLimit: "<<labelRightLimit<<std::endl;
    std::cout<<"label3Limit: "<<label3Limit<<std::endl;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    double sumOfSpectralComponentFour = 0;
    double sumOfSpectralComponentTime = 0;
    std::fstream file;
    file.open("spectralComponent.txt",std::ios::out);
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponentFour,sumOfSpectralComponentTime)
    for(int labelLeft = 0;labelLeft<=labelLeftLimit;labelLeft++){
        if(labelLeft%100==0){
            std::cout<<"labelLeft: "<<labelLeft<<std::endl;
        }
        for(int labelRight = -labelRightLimit;labelRight<=labelRightLimit;labelRight++){
            std::vector<double> emissionPolarAngle = this->calculateEmissionPolarAngle(labelLeft,labelRight);
            //2*emissionPolarAngle
            for(int emissionPolarAngleIndex=0;emissionPolarAngleIndex<emissionPolarAngle.size();emissionPolarAngleIndex++){
                std::complex<double> spectralComponentT =0;
                std::complex<double> spectralComponentX =0;
                std::complex<double> spectralComponentY =0;
                std::complex<double> spectralComponentZ =0;
                for(int label3 = -label3Limit;label3<=label3Limit;label3++){
                    //std::cout<<labelLeft<<","<<labelRight<<","<<label3<<std::endl;
                    std::vector<std::complex<double>> spectralComponent = this->SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngle[emissionPolarAngleIndex]);
                    spectralComponentT += spectralComponent[0];
                    spectralComponentX += spectralComponent[1];
                    spectralComponentY += spectralComponent[2];
                    spectralComponentZ += spectralComponent[3];
                }
                double spectralComponentFour = std::abs(spectralComponentT)*std::abs(spectralComponentT)-std::abs(spectralComponentX)*std::abs(spectralComponentX)-std::abs(spectralComponentY)*std::abs(spectralComponentY)-std::abs(spectralComponentZ)*std::abs(spectralComponentZ);
                double spectralComponentTime = std::abs(spectralComponentT)*std::abs(spectralComponentT);
                file<<"SL: "<<labelLeft<<'\t'<<" SR: "<<labelRight<<'\t'<<"spectralComponentT: "<<std::abs(spectralComponentT)*std::abs(spectralComponentT)<<'\t'<<"spectralComponentX: "<<std::abs(spectralComponentX)*std::abs(spectralComponentX)<<'\t'<<"spectralComponentY: "<<std::abs(spectralComponentY)*std::abs(spectralComponentY)<<'\t'<<"spectralComponentZ: "<<std::abs(spectralComponentZ)*std::abs(spectralComponentZ)<<std::endl;
                //std::cout<<"SL: "<<labelLeft<<'\t'<<" SR: "<<labelRight<<'\t'<<"spectralComponentFour: "<<spectralComponentFour<<'\t'<<"spectralComponentTime: "<<spectralComponentTime<<std::endl;
                sumOfSpectralComponentFour += spectralComponentFour*std::abs(1/(this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*(1/std::tan(emissionPolarAngle[emissionPolarAngleIndex]))-this->getVelocityZPrime()));
                sumOfSpectralComponentTime += spectralComponentTime*std::abs(1/(this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*(1/std::tan(emissionPolarAngle[emissionPolarAngleIndex]))-this->getVelocityZPrime()));
                //std::cout<<"spectral component"<<spectralComponent<<std::endl;
                /*if(std::abs(sumOfSpectralComponentFour)>1){
                    std::cout<<"spectral component"<<sumOfSpectralComponentFour<<std::endl;
                    std::cout<<"labelLeft: "<<labelLeft<<" labelRight: "<<labelRight<<std::endl;
                }*/
            }
        }
    }
    file.close();
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    //print the time with min and sec
    std::cout<<"Time: "<<(time1-time0)/60000000000<<" min"<<std::endl;
    std::cout<<"sumOfSpectralComponentFour: "<<sumOfSpectralComponentFour<<std::endl;
    std::cout<<"sumOfSpectralComponentTime: "<<sumOfSpectralComponentTime<<std::endl;
    double fineStructureConstant = 1.0/137;
    double firstPartOfDifferentialEmissionIntensity = -((fineStructureConstant*electronMass*electronMass*photonEnergy)/(4*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()*residualEnergy))*(this->getEnergy()*this->getEnergy()+residualEnergy*residualEnergy)*sumOfSpectralComponentFour;
    double secondPartOfDifferentialEmissionIntensity = ((fineStructureConstant*electronMass*electronMass*electronMass*electronMass*photonEnergy*photonEnergy*photonEnergy)/(4*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()*this->getEnergy()*this->getEnergy()*residualEnergy))*sumOfSpectralComponentTime;
    this->differentialEmissionIntensity = 2*M_PI*firstPartOfDifferentialEmissionIntensity+secondPartOfDifferentialEmissionIntensity;
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

void BasicRadiationOfElectronInCounterpropagatingLaser::test() {
    std::cout<<"Test:"<<std::endl;
    std::cout<<"=================TEST===================="<<std::endl;
    //std::cout<<"energy/mass"<<this->getEnergy()/electronMass<<std::endl;
    //std::cout<<"parameter1"<<this->getFieldParameter1()<<std::endl;
    //std::cout<<"parameter2"<<this->getFieldParameter2()<<std::endl;
    //std::cout<<"pz/mass"<<this->getInitialMomentumZ()/electronMass<<std::endl;
    //std::vector<double> emissionPolarAngle = this->calculateEmissionPolarAngle(2210,9);
    //std::cout<<"emissionPolarAngle: "<<emissionPolarAngle[0]<<std::endl;
    //std::cout<<"emissionPolarAngle: "<<emissionPolarAngle[1]<<std::endl;
    //double trigonometricCoefficient1 = this->trigonometricCoefficient1(emissionPolarAngle[0]);
    //double trigonometricCoefficient2 = this->trigonometricCoefficient2(emissionPolarAngle[0]);
    //double trigonometricCoefficient3 = this->trigonometricCoefficient3(emissionPolarAngle[0]);
    //double auxiliaryAngle1 = this->auxiliaryAngle1(emissionPolarAngle[0]);
    //double auxiliaryAngle2 = this->auxiliaryAngle2(emissionPolarAngle[0]);
    //std::complex<double> relatedBessel01 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(2216,trigonometricCoefficient1,auxiliaryAngle1);
    //std::complex<double> relatedBessel02 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(2216,trigonometricCoefficient2,auxiliaryAngle2);
    //std::complex<double> relatedBessel03 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(2216,trigonometricCoefficient3,0);
    //std::cout<<"relatedBessel01: "<<relatedBessel01<<std::endl;
    //std::cout<<"relatedBessel02: "<<relatedBessel02<<std::endl;
    //std::cout<<"relatedBessel03: "<<relatedBessel03<<std::endl;
    //std::cout<<BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(0,1,0)<<std::endl;
    std::cout<<"=================TEST===================="<<std::endl;
}

BasicRadiationOfElectronInCounterpropagatingLaser::~BasicRadiationOfElectronInCounterpropagatingLaser() {
    std::cout << "The object has been deleted" << std::endl;
}

