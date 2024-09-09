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
    if(std::abs(argument)<1000&&std::abs(label)<1000){
        return gsl_sf_bessel_Jn(label,argument);
    }
    if(label>=0||((label%2)==0)) {
        return gsl_sf_bessel_Jnu(std::abs(label),argument);
    }else {
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

void BasicRadiationOfElectronInCounterpropagatingLaser::setDifferentialEmissionIntensity(double differentialEmissionIntensity_) {
    this->differentialEmissionIntensity = differentialEmissionIntensity_;
}

std::vector<double> BasicRadiationOfElectronInCounterpropagatingLaser::calculateEmissionPolarAngle(int labelLeft, int labelRight) {
    std::vector<double> emissionPolarAngle={};
    double rho = (1/(this->getEnergy()*energyRatio))*(labelLeft*this->getOmega1()+labelRight*this->getOmega2());
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

std::vector<int> BasicRadiationOfElectronInCounterpropagatingLaser::calculateLabelLimits() {
    int label1Limit = 0;
    int label2Limit = 0;
    int label3Limit = 0;
    double trigonometricCoefficient1Max = this->trigonometricCoefficient1(M_PI/2);
    trigonometricCoefficient1Max = this->trigonometricCoefficient1((std::max(this->getFieldParameter1(),this->getFieldParameter2())+5)*std::sqrt(1-this->getVelocityZPrime()*this->getVelocityZPrime()));
    double trigonometricCoefficient2Max = this->trigonometricCoefficient2(M_PI/2);
    trigonometricCoefficient2Max = this->trigonometricCoefficient2((std::max(this->getFieldParameter1(),this->getFieldParameter2())+5)*std::sqrt(1-this->getVelocityZPrime()*this->getVelocityZPrime()));
    double trigonometricCoefficient3Max = this->trigonometricCoefficient3(0);
    std::cout<<"z1Max: "<<trigonometricCoefficient1Max<<std::endl;
    std::cout<<"z2Max: "<<trigonometricCoefficient2Max<<std::endl;
    std::cout<<"z3Max: "<<trigonometricCoefficient3Max<<std::endl;
    long double timej0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    for(int i = 0;;i++){
        if(std::abs(myBesselFunction(i,trigonometricCoefficient1Max))<1e-25){
            label1Limit = i;
            break;
        }
    }
    for(int i = 0;;i++){
        if(std::abs(myBesselFunction(i,trigonometricCoefficient2Max))<1e-25){
            label2Limit = i;
            break;
        }
    }
    for(int i = 0;;i++){
        if(std::abs(myBesselFunction(i,trigonometricCoefficient3Max))<1e-25){
            label3Limit = i;
            break;
        }
    }
    long double timej1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"Timej: "<<(timej1-timej0)/60000000000<<" min"<<std::endl;
    return {label1Limit,label2Limit,label3Limit};
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
    double auxiliaryAngle1 = this->auxiliaryAngle1(emissionPolarAngle);
    double auxiliaryAngle2 = this->auxiliaryAngle2(emissionPolarAngle);
    double bessel1 = myBesselFunction(label1,trigonometricCoefficient1);
    double bessel2 = myBesselFunction(label2,trigonometricCoefficient2);
    double bessel3 = myBesselFunction(label3,trigonometricCoefficient3);
    double bessel1Plus = myBesselFunction(label1+1,trigonometricCoefficient1);
    double bessel2Plus = myBesselFunction(label2+1,trigonometricCoefficient2);
    double bessel3Plus = myBesselFunction(label3+1,trigonometricCoefficient3);
    double bessel1Minus = myBesselFunction(label1-1,trigonometricCoefficient1);
    double bessel2Minus = myBesselFunction(label2-1,trigonometricCoefficient2);
    double bessel3Minus = myBesselFunction(label3-1,trigonometricCoefficient3);
    //std::cout<<bessel1<<'\t'<<bessel2<<'\t'<<bessel3<<std::endl;
    std::complex<double> zero1 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(bessel1,label1,auxiliaryAngle1);
    std::complex<double> zero2 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(bessel2,label2,auxiliaryAngle2);
    std::complex<double> zero3 = BasicMathFunctionDefinition::relatedBesselFunctionZeroKind(bessel3,label3,0);
    std::complex<double> first1 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(bessel1,bessel1Minus,bessel1Plus,label1,trigonometricCoefficient1,auxiliaryAngle1);
    std::complex<double> first2 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(bessel2,bessel2Minus,bessel2Plus,label2,trigonometricCoefficient2,auxiliaryAngle2);
    std::complex<double> first3 = BasicMathFunctionDefinition::relatedBesselFunctionFirstKind(bessel3,bessel3Minus,bessel3Plus,label3,trigonometricCoefficient3,0);
    std::complex<double> second1 = BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(bessel1,bessel1Minus,bessel1Plus,label1,trigonometricCoefficient1,auxiliaryAngle1);
    std::complex<double> second2 = BasicMathFunctionDefinition::relatedBesselFunctionSecondKind(bessel2,bessel2Minus,bessel2Plus,label2,trigonometricCoefficient2,auxiliaryAngle2);
    std::complex<double> spectralComponentT = (
        zero3*(
            (this->getEnergy()/electronMass)*zero1*zero2+
            ((this->getInitialMomentumX()*omega*this->getFieldParameter1())/(this->getOmega1()*this->getEnergy()))*first1*zero2+
            ((this->getInitialMomentumX()*omega*this->getFieldParameter2())/(this->getOmega2()*this->getEnergy()))*zero1*first2
        )
    );
    std::complex<double> spectralComponentX = (
        zero3*(
            (this->getInitialMomentumX()/electronMass)*zero1*zero2+
            this->getFieldParameter1()*first1*zero2+
            this->getFieldParameter2()*zero1*first2
        )
    );
    std::complex<double> spectralComponentY = (zero3*(
            this->getFieldParameter1()*second1*zero2+
            this->getFieldParameter2()*zero1*second2
        )
    );
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
    std::complex<double> spectralComponentBase = zero1*zero2*zero3;
    std::vector<std::complex<double>> spectralComponent = {spectralComponentT,spectralComponentX,spectralComponentY,spectralComponentZ,spectralComponentBase};
    return spectralComponent;
}

void BasicRadiationOfElectronInCounterpropagatingLaser::calculateDifferentialEmissionIntensity() {
    std::cout<<"this is the function in BasicRadiationOfElectronInCounterpropagatingLaser"<<std::endl;
    std::vector<int> labelLimits = this->calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = 20000;
    //labelRightLimit = 30;
    label3Limit = 10;
    std::cout<<"labelLeftLimit: "<<labelLeftLimit<<std::endl;
    std::cout<<"labelRightLimit: "<<labelRightLimit<<std::endl;
    std::cout<<"label3Limit: "<<label3Limit<<std::endl;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    double sumOfSpectralComponentFour = 0;
    double sumOfSpectralComponentTime = 0;
    //std::fstream file;
    //file.open("spectralComponent1.txt",std::ios::out);
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
                    //file<<labelLeft<<","<<labelRight<<","<<label3<<std::endl;
                    std::vector<std::complex<double>> spectralComponent = this->SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngle[emissionPolarAngleIndex]);
                    //file<<spectralComponent[0]<<","<<spectralComponent[1]<<","<<spectralComponent[2]<<","<<spectralComponent[3]<<spectralComponent[4]<<std::endl;
                    spectralComponentT += spectralComponent[0];
                    spectralComponentX += spectralComponent[1];
                    spectralComponentY += spectralComponent[2];
                    spectralComponentZ += spectralComponent[3];
                }
                double spectralComponentFour = std::abs(spectralComponentT)*std::abs(spectralComponentT)-std::abs(spectralComponentX)*std::abs(spectralComponentX)-std::abs(spectralComponentY)*std::abs(spectralComponentY)-std::abs(spectralComponentZ)*std::abs(spectralComponentZ);
                double spectralComponentTime = std::abs(spectralComponentT)*std::abs(spectralComponentT);
                //file<<"SL: "<<labelLeft<<'\t'<<" SR: "<<labelRight<<'\t'<<"spectralComponentT: "<<std::abs(spectralComponentT)*std::abs(spectralComponentT)<<'\t'<<"spectralComponentX: "<<std::abs(spectralComponentX)*std::abs(spectralComponentX)<<'\t'<<"spectralComponentY: "<<std::abs(spectralComponentY)*std::abs(spectralComponentY)<<'\t'<<"spectralComponentZ: "<<std::abs(spectralComponentZ)*std::abs(spectralComponentZ)<<std::endl;
                //std::cout<<"SL: "<<labelLeft<<'\t'<<" SR: "<<labelRight<<'\t'<<"spectralComponentFour: "<<spectralComponentFour<<'\t'<<"spectralComponentTime: "<<spectralComponentTime<<std::endl;
                sumOfSpectralComponentFour += spectralComponentFour*std::abs(1/(this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*(1/std::tan(emissionPolarAngle[emissionPolarAngleIndex]))-this->getVelocityZPrime()));
                sumOfSpectralComponentTime += spectralComponentTime*std::abs(1/(this->getVelocityXPrime()*std::cos(this->emissionAzimuthalAngle)*(1/std::tan(emissionPolarAngle[emissionPolarAngleIndex]))-this->getVelocityZPrime()));
            }
        }
    }
    //file.close();
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    //print the time with min and sec
    std::cout<<"Time: "<<(time1-time0)/60000000000<<" min"<<std::endl;
    std::cout<<"sumOfSpectralComponentFour: "<<sumOfSpectralComponentFour<<std::endl;
    std::cout<<"sumOfSpectralComponentTime: "<<sumOfSpectralComponentTime<<std::endl;
    double fineStructureConstant = 1.0/137;
    double firstPartOfDifferentialEmissionIntensity = -((fineStructureConstant*electronMass*electronMass*photonEnergy)/(4*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()*residualEnergy))*(this->getEnergy()*this->getEnergy()+residualEnergy*residualEnergy)*sumOfSpectralComponentFour;
    double secondPartOfDifferentialEmissionIntensity = ((fineStructureConstant*electronMass*electronMass*electronMass*electronMass*photonEnergy*photonEnergy*photonEnergy)/(4*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()*this->getEnergy()*this->getEnergy()*residualEnergy))*sumOfSpectralComponentTime;
    this->differentialEmissionIntensity = firstPartOfDifferentialEmissionIntensity+secondPartOfDifferentialEmissionIntensity;
    std::cout<<this->differentialEmissionIntensity<<std::endl;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getPhotonEnergy() {
    return this->photonEnergy;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getEmissionAzimuthalAngle() {
    return this->emissionAzimuthalAngle;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getDifferentialEmissionIntensity() {
    return this->differentialEmissionIntensity;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getResidualEnergy() {
    return this->residualEnergy;
}

double BasicRadiationOfElectronInCounterpropagatingLaser::getEnergyRatio() {
    return this->energyRatio;
}

void BasicRadiationOfElectronInCounterpropagatingLaser::test() {
    std::cout<<"Test:"<<std::endl;
    std::cout<<"=================TEST===================="<<std::endl;
    std::cout<<"=================TEST===================="<<std::endl;
}

BasicRadiationOfElectronInCounterpropagatingLaser::~BasicRadiationOfElectronInCounterpropagatingLaser() {
    std::cout << "The object has been deleted" << std::endl;
}

