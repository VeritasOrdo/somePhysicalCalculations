#include "RadiationWithSpinAndPolarzation.h"
#include <chrono>
#include <omp.h>
#include <fstream>  

template <typename T>
T sinc_unnormalized(T x) {
    // 当 x 的绝对值非常接近0时，返回1.0
    // 使用 epsilon 比直接用 x == 0.0 更安全，可以处理极小的非零值
    if (std::abs(x) < std::numeric_limits<T>::epsilon()) {
        return T(1.0);
    }
    return std::sin(x) / x;
}

RadiationWithSpinAndPolarzation::RadiationWithSpinAndPolarzation(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime,double photonEnergy,double emissionAzimuthalAngle,double spinIncident,double spinEmission,double polarizationAlpha,double polarizationBeta,double axisOfIncidentAzimuthalAngleOfElectronSpin,double axisOfIncidentPolarAngleOfElectronSpin,double axisOfEmissionAzimuthalAngleOfElectronSpin,double axisOfEmissionPolarAngleOfElectronSpin,double rotationDirection1, double rotationDirection2,double emissionPolarAngleMin,double emissionPolarAngleMax) : BasicRadiationOfElectronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,properTime,photonEnergy,emissionAzimuthalAngle,rotationDirection1,rotationDirection2) {
    this->spinIncident = spinIncident;
    this->spinEmission = spinEmission;
    this->polarizationAlpha = polarizationAlpha;
    this->polarizationBeta = polarizationBeta;
    this->axisOfIncidentAzimuthalAngleOfElectronSpin = axisOfIncidentAzimuthalAngleOfElectronSpin;
    this->axisOfIncidentPolarAngleOfElectronSpin = axisOfIncidentPolarAngleOfElectronSpin;
    this->axisOfEmissionAzimuthalAngleOfElectronSpin = axisOfEmissionAzimuthalAngleOfElectronSpin;
    this->axisOfEmissionPolarAngleOfElectronSpin = axisOfEmissionPolarAngleOfElectronSpin;
    this->incidentOrientationAxis = Dimension3Vector<double>(std::cos(axisOfIncidentAzimuthalAngleOfElectronSpin)*std::sin(axisOfIncidentPolarAngleOfElectronSpin),std::sin(axisOfIncidentAzimuthalAngleOfElectronSpin)*std::sin(axisOfIncidentPolarAngleOfElectronSpin),std::cos(axisOfIncidentPolarAngleOfElectronSpin));
    this->emissionOrientationAxis = Dimension3Vector<double>(std::cos(axisOfEmissionAzimuthalAngleOfElectronSpin)*std::sin(axisOfEmissionPolarAngleOfElectronSpin),std::sin(axisOfEmissionAzimuthalAngleOfElectronSpin)*std::sin(axisOfEmissionPolarAngleOfElectronSpin),std::cos(axisOfEmissionPolarAngleOfElectronSpin));
    this->combinedIncidentOrientationAxis = incidentOrientationAxis*this->spinIncident*2.0;
    this->combinedEmissionOrientationAxis = emissionOrientationAxis*this->spinEmission*2.0;
    this->emissionPolarAngleMax = emissionPolarAngleMax;
    this->emissionPolarAngleMin = emissionPolarAngleMin;
    this->stokesParameter = {0,0,0,0};
}

void RadiationWithSpinAndPolarzation::calculateDifferentialEmissionIntensity() {
    std::cout<<"this is the function in RadiationWithSpinAndPolarzation"<<std::endl;
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = 10;
    //labelRightLimit = 30;
    //label3Limit = 10;
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    //std::fstream file;
    //file.open("spectralComponent.txt", std::ios::out);
    double sumOfSpectralComponent = 0;
    //double sumOfSpectralComponentImag = 0;
    //double sumOfSpectralComponentTest = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"max thread: "<<omp_get_max_threads()<<std::endl;
    std::cout << "label left limit min: " << -std::min(std::max(labelLeftLimit/100,500),40000) << std::endl;
    int count = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponent) reduction(+:count)
    for(int labelLeft = -std::min(std::max(labelLeftLimit/100,500),40000); labelLeft <=labelLeftLimit; labelLeft++) {
        if(labelLeft%100==0){
            std::cout<<"SL: "<<labelLeft<<std::endl;
        }
        for(int labelRight = -labelRightLimit; labelRight <=labelRightLimit; labelRight++) {
            /*if(labelRight%100==0){
                std::cout<<"SR: "<<labelRight<<std::endl;
            }*/
            std::vector<double> emissionPolarAngles = calculateEmissionPolarAngle(labelLeft,labelRight);
            std::complex<double> sumOfComponentA=0;
            std::complex<double> sumOfComponentAX=0;
            std::complex<double> sumOfComponentAY=0;
            Dimension3Vector<std::complex<double>> sumOfComponentB = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBX = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBY = Dimension3Vector<std::complex<double>>(0,0,0);
            for(int label3 = -label3Limit; label3 <=label3Limit; label3++) {
                for(int labelOfEmissionPolarAngles = 0; labelOfEmissionPolarAngles < emissionPolarAngles.size(); labelOfEmissionPolarAngles++) {
                    if((emissionPolarAngles[labelOfEmissionPolarAngles]<emissionPolarAngleMin)||(emissionPolarAngles[labelOfEmissionPolarAngles]>emissionPolarAngleMax)){
                        continue;
                    }
                    count++;
                    std::vector<std::complex<double>> spectralComponent = SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngles[labelOfEmissionPolarAngles]);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1],spectralComponent[2],spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    //Dimension3Vector<std::complex<double>> spectralComponentPhotonEmissionVector = ((photonEmissionVector*spectralComponent[4])/this->electronMass)*this->getPhotonEnergy();
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVector = polarizationVectorBase1*std::cos(polarizationAlpha)+polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    //polarizationVector = ((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized()+(((photonEmissionVector)^((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized())).normalized())*std::complex<double>(0,1))*(1.0/sqrt(2.0));
                    Dimension3Vector<std::complex<double>> polarizationVectorX = polarizationVectorBase1*std::cos(polarizationAlpha);
                    Dimension3Vector<std::complex<double>> polarizationVectorY = polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    Dimension3Vector<std::complex<double>> polarizationVectorConjugate = polarizationVector.conjugate();
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy()+this->electronMass)/(this->getEnergy()+this->electronMass));
                    double secondRatioOfEnergy = 1/firstRatioOfEnergy;
                    //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<" S3: "<<label3<<std::endl;
                    sumOfComponentA += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorConjugate*spectralComponent3D);
                    //sumOfComponentB += ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorConjugate^(spectralComponent3D-spectralComponentPhotonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    sumOfComponentB += ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-((polarizationVectorConjugate^(spectralComponent3D*this->getEnergy()-photonEmissionVector*this->getPhotonEnergy()*spectralComponent[0])*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))));
                }
            }
            std::complex<double> componentASpinCoefficient = (
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))                    )
            );
            Dimension3Vector<std::complex<double>> componentBSpinCoefficient = Dimension3Vector<std::complex<double>>(
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,1)
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,-1)
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )-
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )
            );
            //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<std::endl;
            //file<<componentASpinCoefficient<<'\t'<<componentBSpinCoefficient<<std::endl;
            //file<<sumOfComponentA*componentASpinCoefficient<<'\t'<<sumOfComponentB*componentBSpinCoefficient<<std::endl;
            //std::complex<double> sumOfSpectralComponentAmplitude = sumOfComponentA*componentASpinCoefficient+sumOfComponentB*componentBSpinCoefficient*std::complex<double>(0,1);
            //sumOfSpectralComponentTest += std::abs(sumOfSpectralComponentAmplitude)*std::abs(sumOfSpectralComponentAmplitude);
            sumOfSpectralComponent += (
                (0.5+(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(std::conj(sumOfComponentA)*sumOfComponentA)+
                (0.5-(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(sumOfComponentB.conjugate()*sumOfComponentB)+
                (sumOfComponentB*combinedEmissionOrientationAxis)*(sumOfComponentB.conjugate()*combinedIncidentOrientationAxis)-
                (-(sumOfComponentB.conjugate()^sumOfComponentB)+(sumOfComponentB.conjugate()*sumOfComponentA)+(sumOfComponentB*std::conj(sumOfComponentA)))*(combinedEmissionOrientationAxis^combinedIncidentOrientationAxis)*(0.5)+
                (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedEmissionOrientationAxis-combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)+
                ((sumOfComponentB*(std::conj(sumOfComponentA)))-(sumOfComponentB.conjugate()*sumOfComponentA))*(combinedEmissionOrientationAxis+combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)
            ).real();
        }
    }
    std::cout<<"count: "<<count<<std::endl;
    // theta+5/gamma
    //file.close();
    std::cout<<"sumOfSpectralComponentReal: "<<sumOfSpectralComponent<<std::endl;
    //std::cout<<"sumOfSpectralComponentImag: "<<sumOfSpectralComponentImag<<std::endl;
    //std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    double fineStructureConstant = 1.0/137;
    this->setDifferentialEmissionIntensity(((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponent);
    std::cout<<"differentialEmissionIntensity: "<<this->getDifferentialEmissionIntensity()<<std::endl;
    //std::cout<<"differentialEmissionIntensityTest: "<<((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponentTest<<std::endl;
}

void RadiationWithSpinAndPolarzation::calculateDifferentialEmissionIntensityWithDoubledLabel() {
    std::cout<<"this is the function in RadiationWithSpinAndPolarzation"<<std::endl;
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = 10;
    //labelRightLimit = 30;
    //label3Limit = 10;
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    //std::fstream file;
    //file.open("spectralComponent.txt", std::ios::out);
    double sumOfSpectralComponent = 0;
    //double sumOfSpectralComponentImag = 0;
    //double sumOfSpectralComponentTest = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"max thread: "<<omp_get_max_threads()<<std::endl;
    std::cout << "label left limit min: " << -std::min(std::max(labelLeftLimit/100,500),40000) << std::endl;
    int count = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponent) reduction(+:count)
    for(int labelLeft = -2*std::min(std::max(labelLeftLimit/100,500),40000); labelLeft <=2*labelLeftLimit; labelLeft++) {
        if(labelLeft%100==0){
            std::cout<<"SL: "<<labelLeft<<std::endl;
        }
        for(int labelRight = -2*labelRightLimit; labelRight <=2*labelRightLimit; labelRight++) {
            /*if(labelRight%100==0){
                std::cout<<"SR: "<<labelRight<<std::endl;
            }*/
            std::vector<double> emissionPolarAngles = calculateEmissionPolarAngle(labelLeft,labelRight);
            std::complex<double> sumOfComponentA=0;
            std::complex<double> sumOfComponentAX=0;
            std::complex<double> sumOfComponentAY=0;
            Dimension3Vector<std::complex<double>> sumOfComponentB = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBX = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBY = Dimension3Vector<std::complex<double>>(0,0,0);
            for(int label3 = -2*label3Limit; label3 <=2*label3Limit; label3++) {
                for(int labelOfEmissionPolarAngles = 0; labelOfEmissionPolarAngles < emissionPolarAngles.size(); labelOfEmissionPolarAngles++) {
                    if((emissionPolarAngles[labelOfEmissionPolarAngles]<emissionPolarAngleMin)||(emissionPolarAngles[labelOfEmissionPolarAngles]>emissionPolarAngleMax)){
                        continue;
                    }
                    count++;
                    std::vector<std::complex<double>> spectralComponent = SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngles[labelOfEmissionPolarAngles]);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1],spectralComponent[2],spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    //Dimension3Vector<std::complex<double>> spectralComponentPhotonEmissionVector = ((photonEmissionVector*spectralComponent[4])/this->electronMass)*this->getPhotonEnergy();
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVector = polarizationVectorBase1*std::cos(polarizationAlpha)+polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    //polarizationVector = ((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized()+(((photonEmissionVector)^((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized())).normalized())*std::complex<double>(0,1))*(1.0/sqrt(2.0));
                    Dimension3Vector<std::complex<double>> polarizationVectorX = polarizationVectorBase1*std::cos(polarizationAlpha);
                    Dimension3Vector<std::complex<double>> polarizationVectorY = polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    Dimension3Vector<std::complex<double>> polarizationVectorConjugate = polarizationVector.conjugate();
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy()+this->electronMass)/(this->getEnergy()+this->electronMass));
                    double secondRatioOfEnergy = 1/firstRatioOfEnergy;
                    //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<" S3: "<<label3<<std::endl;
                    sumOfComponentA += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorConjugate*spectralComponent3D);
                    //sumOfComponentB += ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorConjugate^(spectralComponent3D-spectralComponentPhotonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    sumOfComponentB += ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-((polarizationVectorConjugate^(spectralComponent3D*this->getEnergy()-photonEmissionVector*this->getPhotonEnergy()*spectralComponent[0])*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))));
                }
            }
            std::complex<double> componentASpinCoefficient = (
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))                    )
            );
            Dimension3Vector<std::complex<double>> componentBSpinCoefficient = Dimension3Vector<std::complex<double>>(
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,1)
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,-1)
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )-
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )
            );
            sumOfSpectralComponent += (
                (0.5+(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(std::conj(sumOfComponentA)*sumOfComponentA)+
                (0.5-(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(sumOfComponentB.conjugate()*sumOfComponentB)+
                (sumOfComponentB*combinedEmissionOrientationAxis)*(sumOfComponentB.conjugate()*combinedIncidentOrientationAxis)-
                (-(sumOfComponentB.conjugate()^sumOfComponentB)+(sumOfComponentB.conjugate()*sumOfComponentA)+(sumOfComponentB*std::conj(sumOfComponentA)))*(combinedEmissionOrientationAxis^combinedIncidentOrientationAxis)*(0.5)+
                (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedEmissionOrientationAxis-combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)+
                ((sumOfComponentB*(std::conj(sumOfComponentA)))-(sumOfComponentB.conjugate()*sumOfComponentA))*(combinedEmissionOrientationAxis+combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)
            ).real();
        }
    }
    std::cout<<"count: "<<count<<std::endl;
    // theta+5/gamma
    //file.close();
    std::cout<<"sumOfSpectralComponentReal: "<<sumOfSpectralComponent<<std::endl;
    //std::cout<<"sumOfSpectralComponentImag: "<<sumOfSpectralComponentImag<<std::endl;
    //std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    double fineStructureConstant = 1.0/137;
    this->setDifferentialEmissionIntensity(((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponent);
    std::cout<<"differentialEmissionIntensity: "<<this->getDifferentialEmissionIntensity()<<std::endl;
    //std::cout<<"differentialEmissionIntensityTest: "<<((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponentTest<<std::endl;
}

void RadiationWithSpinAndPolarzation::calculateVortexDifferentialEmissionIntensity(double angularQuantumNumber, double polarizationParameter, size_t azimuthalAngleDivisions, double emissionPolarAngle, double nForSinc){
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = labelLeftLimit*2;
    //labelRightLimit = labelRightLimit*2;
    //label3Limit = label3Limit*2;
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    double sumOfSpectralComponentReal = 0;
    double sumOfSpectralComponentImag = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"max thread: "<<omp_get_max_threads()<<std::endl;
    std::cout << "label left limit min: " << -std::min(std::max(labelLeftLimit/100,500),40000) << std::endl;
    double azimuthalAngleStep = 2*M_PI/azimuthalAngleDivisions;
    double angleRelatedCoeffcient = ((this->getPhotonEnergy()/(this->getEnergy()-this->getPhotonEnergy()))*this->getEnergy()*this->getEnergy()/electronMass)*(1-this->getVelocityZPrime()*std::cos(emissionPolarAngle));
    
    // 在临界区外进行自旋系数计算
    std::complex<double> componentASpinCoefficient = (
        (
            (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
        )+
        (
            (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
        )
    );
            
    Dimension3Vector<std::complex<double>> componentBSpinCoefficient = Dimension3Vector<std::complex<double>>(
        (
            (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
        )+
        (
            (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
        ),
        (
            (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
            std::complex<double>(0,1)
        )+
        (
            (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
            std::complex<double>(0,-1)
        ),
        (
            (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
        )-
        (
            (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
            (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
        )
    );

    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponentReal) reduction(+:sumOfSpectralComponentImag)
    for (int azimuthalAngleLabel = 0; azimuthalAngleLabel < azimuthalAngleDivisions; azimuthalAngleLabel++)
    {
        double azimuthalAngle = azimuthalAngleLabel * azimuthalAngleStep;

        RadiationWithSpinAndPolarzation localRadiationWithSpinAndPolarzation(*this);
        localRadiationWithSpinAndPolarzation.setEmissionAzimuthalAngle(azimuthalAngle);

        if (azimuthalAngleLabel % 10 == 0)
        {
            #pragma omp critical
            {
                std::cout << "Azimuthal Angle Label: " << azimuthalAngleLabel << std::endl;
            }
        }

        for (int labelLeft = -std::min(std::max(labelLeftLimit / 100, 500), 40000); labelLeft <= labelLeftLimit; labelLeft++)
        {

            int labelRightMidium = int((angleRelatedCoeffcient / (this->getEnergy() / electronMass) - (labelLeft * this->getOmega1())) / this->getOmega2());

            for (int labelRight = labelRightMidium - labelRightLimit; labelRight <= labelRightMidium + labelRightLimit; labelRight++)
            {
                std::complex<double> integralOfSpectralComponent = 0;

                double labelRelatedCoeffcient = -(this->getEnergy() / electronMass) * (labelLeft * this->getOmega1() + labelRight * this->getOmega2());
                double deltaCoeffcient = labelRelatedCoeffcient + angleRelatedCoeffcient;
                if(labelRelatedCoeffcient + angleRelatedCoeffcient > nForSinc){
                    continue;
                }
                double deltaReplacedSinc = sinc_unnormalized(deltaCoeffcient)*sinc_unnormalized(deltaCoeffcient/nForSinc);

                std::complex<double> sumOfComponentA = 0;
                Dimension3Vector<std::complex<double>> sumOfComponentB = Dimension3Vector<std::complex<double>>(0, 0, 0);

                for (int label3 = -label3Limit; label3 <= label3Limit; label3++)
                {
                    std::vector<std::complex<double>> spectralComponent = localRadiationWithSpinAndPolarzation.SpectralComponent(labelLeft, labelRight, label3, emissionPolarAngle);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1], spectralComponent[2], spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngle) * std::cos(azimuthalAngle), std::sin(emissionPolarAngle) * std::sin(azimuthalAngle), std::cos(emissionPolarAngle));
                    std::complex<double> I = std::complex<double>(0, 1);
                    Dimension3Vector<std::complex<double>> vortexBasePlus = Dimension3Vector<std::complex<double>>(1, polarizationParameter * I, 0) * (-polarizationParameter / std::sqrt(2));
                    Dimension3Vector<std::complex<double>> vortexBaseMinus = Dimension3Vector<std::complex<double>>(1, -polarizationParameter * I, 0) * (polarizationParameter / std::sqrt(2));
                    Dimension3Vector<std::complex<double>> vortexBaseZ(0, 0, 1);
                    Dimension3Vector<std::complex<double>> polarizationVector = ((vortexBaseMinus * std::exp(I * polarizationParameter * azimuthalAngle) * std::pow(std::sin(emissionPolarAngle / 2), 2)) +
                                                                                 (vortexBasePlus * std::exp(-I * polarizationParameter * azimuthalAngle) * std::pow(std::cos(emissionPolarAngle / 2), 2)) +
                                                                                 (vortexBaseZ * std::sin(emissionPolarAngle) * (polarizationParameter / 2)));
                    Dimension3Vector<std::complex<double>> polarizationVectorConjugate = polarizationVector.conjugate();
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy() + this->electronMass) / (this->getEnergy() + this->electronMass));
                    double secondRatioOfEnergy = 1 / firstRatioOfEnergy;

                    sumOfComponentA += ((this->getEnergy()) / (2 * std::sqrt(this->getEnergy() * this->getResidualEnergy()))) * (firstRatioOfEnergy + secondRatioOfEnergy) * (polarizationVectorConjugate * spectralComponent3D);
                    sumOfComponentB += ((polarizationVectorConjugate ^ spectralComponent3D) * firstRatioOfEnergy * (this->getEnergy()) / (2 * std::sqrt(this->getEnergy() * this->getResidualEnergy()))) - ((polarizationVectorConjugate ^ (spectralComponent3D * this->getEnergy() - photonEmissionVector * this->getPhotonEnergy() * spectralComponent[0]) * secondRatioOfEnergy * (1.0 / (2 * std::sqrt(this->getEnergy() * this->getResidualEnergy())))));
                }

                std::complex<double> spectralComponentAmplitude = sumOfComponentA * componentASpinCoefficient + sumOfComponentB * componentBSpinCoefficient * std::complex<double>(0, 1);
                integralOfSpectralComponent = spectralComponentAmplitude * std::pow(std::complex<double>(0, 1), angularQuantumNumber) * std::exp(std::complex<double>(0, -1) * angularQuantumNumber * azimuthalAngle) * (1.0 / std::pow(std::sqrt(2 * M_PI), 3)) * std::sqrt(this->getPhotonEnergy() * std::sin(emissionPolarAngle)) * azimuthalAngleStep;

                sumOfSpectralComponentReal += (integralOfSpectralComponent * deltaReplacedSinc).real();
                sumOfSpectralComponentImag += (integralOfSpectralComponent * deltaReplacedSinc).imag();
            }
        }
    }

    //std::cout<<"sumOfSpectralComponentReal: "<<sumOfSpectralComponent<<std::endl;
    //std::cout<<"sumOfSpectralComponentImag: "<<sumOfSpectralComponentImag<<std::endl;
    //std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    double fineStructureConstant = 1.0/137;
    std::complex<double> sumOfSpectralComponent = std::complex<double>(sumOfSpectralComponentReal, sumOfSpectralComponentImag);
    this->setDifferentialEmissionIntensity(std::norm(sumOfSpectralComponent));
    std::cout<<"differentialEmissionIntensity: "<<this->getDifferentialEmissionIntensity()<<std::endl;
}

void RadiationWithSpinAndPolarzation::calculateStokesParameter() {
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = 10;
    //labelRightLimit = 30;
    //label3Limit = 10;
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    //std::fstream file;
    //file.open("spectralComponent.txt", std::ios::out);
    double sumOfSpectralComponentStokesParameterI = 0;
    double sumOfSpectralComponentStokesParameterQ = 0;
    double sumOfSpectralComponentStokesParameterU = 0;
    double sumOfSpectralComponentStokesParameterV = 0;
    //double sumOfSpectralComponentImag = 0;
    //double sumOfSpectralComponentTest = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"max thread: "<<omp_get_max_threads()<<std::endl;
    std::cout << "label left limit min: " << -std::min(std::max(labelLeftLimit/100,500),40000) << std::endl;
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponentStokesParameterI) reduction(+:sumOfSpectralComponentStokesParameterQ) reduction(+:sumOfSpectralComponentStokesParameterU) reduction(+:sumOfSpectralComponentStokesParameterV)
    for(int labelLeft = -std::min(std::max(labelLeftLimit/100,500),40000); labelLeft <=labelLeftLimit; labelLeft++) {
        /*if(labelLeft%100==0){
            std::cout<<"SL: "<<labelLeft<<std::endl;
        }*/
        for(int labelRight = -labelRightLimit; labelRight <=labelRightLimit; labelRight++) {
            /*if(labelRight%100==0){
                std::cout<<"SR: "<<labelRight<<std::endl;
            }*/
            std::vector<double> emissionPolarAngles = calculateEmissionPolarAngle(labelLeft,labelRight);
            std::complex<double> sumOfComponentAX=0;
            std::complex<double> sumOfComponentAY=0;
            Dimension3Vector<std::complex<double>> sumOfComponentBX = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBY = Dimension3Vector<std::complex<double>>(0,0,0);
            for(int label3 = -label3Limit; label3 <=label3Limit; label3++) {
                for(int labelOfEmissionPolarAngles = 0; labelOfEmissionPolarAngles < emissionPolarAngles.size(); labelOfEmissionPolarAngles++) {
                    if((emissionPolarAngles[labelOfEmissionPolarAngles]<emissionPolarAngleMin)||(emissionPolarAngles[labelOfEmissionPolarAngles]>emissionPolarAngleMax)){
                        continue;
                    }
                    std::vector<std::complex<double>> spectralComponent = SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngles[labelOfEmissionPolarAngles]);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1],spectralComponent[2],spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    //Dimension3Vector<std::complex<double>> spectralComponentPhoton
                    //EmissionVector = ((photonEmissionVector*spectralComponent[4])/this->electronMass)*this->getPhotonEnergy();
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVectorX = polarizationVectorBase1;
                    Dimension3Vector<std::complex<double>> polarizationVectorY = polarizationVectorBase2;
                    Dimension3Vector<std::complex<double>> polarizationVectorXConjugate = polarizationVectorX.conjugate();
                    Dimension3Vector<std::complex<double>> polarizationVectorYConjugate = polarizationVectorY.conjugate();
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy()+this->electronMass)/(this->getEnergy()+this->electronMass));
                    double secondRatioOfEnergy = 1/firstRatioOfEnergy;
                    //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<" S3: "<<label3<<std::endl;
                    sumOfComponentAX += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorXConjugate*spectralComponent3D);
                    sumOfComponentAY += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorYConjugate*spectralComponent3D);
                    sumOfComponentBX += ((polarizationVectorXConjugate^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-((polarizationVectorXConjugate^(spectralComponent3D*this->getEnergy()-photonEmissionVector*this->getPhotonEnergy()*spectralComponent[0])*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))));
                    sumOfComponentBY += ((polarizationVectorYConjugate^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-((polarizationVectorYConjugate^(spectralComponent3D*this->getEnergy()-photonEmissionVector*this->getPhotonEnergy()*spectralComponent[0])*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))));
                }
            }
            std::complex<double> componentASpinCoefficient = (
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))                    )
            );
            Dimension3Vector<std::complex<double>> componentBSpinCoefficient = Dimension3Vector<std::complex<double>>(
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,1)
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,-1)
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )-
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )
            );
            //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<std::endl;
            //file<<componentASpinCoefficient<<'\t'<<componentBSpinCoefficient<<std::endl;
            //file<<sumOfComponentA*componentASpinCoefficient<<'\t'<<sumOfComponentB*componentBSpinCoefficient<<std::endl;
            std::complex<double> sumOfSpectralComponentAmplitudeX = sumOfComponentAX*componentASpinCoefficient+sumOfComponentBX*componentBSpinCoefficient*std::complex<double>(0,1);
            std::complex<double> sumOfSpectralComponentAmplitudeY = sumOfComponentAY*componentASpinCoefficient+sumOfComponentBY*componentBSpinCoefficient*std::complex<double>(0,1);
            /*if(sumOfSpectralComponentAmplitudeX.real()!=0){
                std::cout<<"sumOfSpectralComponentAmplitudeX: "<<sumOfSpectralComponentAmplitudeX<<std::endl;
                std::cout<<"sumOfSpectralComponentAmplitudeY: "<<sumOfSpectralComponentAmplitudeY<<std::endl;
                std::cout<<"sumOfComponentAY: "<<sumOfComponentAY<<std::endl;
                std::cout<<"sumOfComponentBY: "<<sumOfComponentBY<<std::endl;
            }*/
            //sumOfSpectralComponentTest += std::abs(sumOfSpectralComponentAmplitude)*std::abs(sumOfSpectralComponentAmplitude);
            sumOfSpectralComponentStokesParameterI += std::abs(sumOfSpectralComponentAmplitudeX)*std::abs(sumOfSpectralComponentAmplitudeX)+std::abs(sumOfSpectralComponentAmplitudeY)*std::abs(sumOfSpectralComponentAmplitudeY);
            sumOfSpectralComponentStokesParameterQ += std::abs(sumOfSpectralComponentAmplitudeX)*std::abs(sumOfSpectralComponentAmplitudeX)-std::abs(sumOfSpectralComponentAmplitudeY)*std::abs(sumOfSpectralComponentAmplitudeY);
            sumOfSpectralComponentStokesParameterU += 2.0*std::real(sumOfSpectralComponentAmplitudeX*std::conj(sumOfSpectralComponentAmplitudeY));
            sumOfSpectralComponentStokesParameterV += -2.0*std::imag(sumOfSpectralComponentAmplitudeX*std::conj(sumOfSpectralComponentAmplitudeY));
        }
    }
    // theta+5/gamma
    //file.close();
    std::cout<<"sumOfSpectralComponentStokesParameterI: "<<sumOfSpectralComponentStokesParameterI<<std::endl;
    std::cout<<"sumOfSpectralComponentStokesParameterQ: "<<sumOfSpectralComponentStokesParameterQ<<std::endl;
    std::cout<<"sumOfSpectralComponentStokesParameterU: "<<sumOfSpectralComponentStokesParameterU<<std::endl;
    std::cout<<"sumOfSpectralComponentStokesParameterV: "<<sumOfSpectralComponentStokesParameterV<<std::endl;
    //std::cout<<"sumOfSpectralComponentImag: "<<sumOfSpectralComponentImag<<std::endl;
    //std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    stokesParameter = {sumOfSpectralComponentStokesParameterI,sumOfSpectralComponentStokesParameterQ,sumOfSpectralComponentStokesParameterU,sumOfSpectralComponentStokesParameterV};
    //std::cout<<"differentialEmissionIntensity: "<<this->getDifferentialEmissionIntensity()<<std::endl;
    //std::cout<<"differentialEmissionIntensityTest: "<<((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponentTest<<std::endl;
}

std::vector<double> RadiationWithSpinAndPolarzation::calculateSixTermsOfDifferentialEmissionIntensity() {
    double sixTermsOfDifferentialEmissionIntensity0 = 0;
    double sixTermsOfDifferentialEmissionIntensity1 = 0;
    double sixTermsOfDifferentialEmissionIntensity2 = 0;
    double sixTermsOfDifferentialEmissionIntensity3 = 0;
    double sixTermsOfDifferentialEmissionIntensity4 = 0;
    double sixTermsOfDifferentialEmissionIntensity5 = 0;
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = 10;
    //labelRightLimit = 30;
    //label3Limit = 10;
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout<<"max thread: "<<omp_get_max_threads()<<std::endl;
    std::cout << "label left limit min: " << -std::min(std::max(labelLeftLimit/100,500),40000) << std::endl;
    #pragma omp parallel for schedule(dynamic) reduction(+:sixTermsOfDifferentialEmissionIntensity0) reduction(+:sixTermsOfDifferentialEmissionIntensity1) reduction(+:sixTermsOfDifferentialEmissionIntensity2) reduction(+:sixTermsOfDifferentialEmissionIntensity3) reduction(+:sixTermsOfDifferentialEmissionIntensity4) reduction(+:sixTermsOfDifferentialEmissionIntensity5)
    for(int labelLeft = -std::min(std::max(labelLeftLimit/100,500),40000); labelLeft <=labelLeftLimit; labelLeft++) {
        /*if(labelLeft%100==0){
            std::cout<<"SL: "<<labelLeft<<std::endl;
        }*/
        for(int labelRight = -labelRightLimit; labelRight <=labelRightLimit; labelRight++) {
            /*if(labelRight%100==0){
                std::cout<<"SR: "<<labelRight<<std::endl;
            }*/
            std::vector<double> emissionPolarAngles = calculateEmissionPolarAngle(labelLeft,labelRight);
            std::complex<double> sumOfComponentA=0;
            Dimension3Vector<std::complex<double>> sumOfComponentB = Dimension3Vector<std::complex<double>>(0,0,0);
            for(int label3 = -label3Limit; label3 <=label3Limit; label3++) {
                for(int labelOfEmissionPolarAngles = 0; labelOfEmissionPolarAngles < emissionPolarAngles.size(); labelOfEmissionPolarAngles++) {
                    std::vector<std::complex<double>> spectralComponent = SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngles[labelOfEmissionPolarAngles]);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1],spectralComponent[2],spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    //Dimension3Vector<std::complex<double>> spectralComponentPhoton
                    //EmissionVector = ((photonEmissionVector*spectralComponent[4])/this->electronMass)*this->getPhotonEnergy();
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVector = polarizationVectorBase1*std::cos(polarizationAlpha)+polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy()+this->electronMass)/(this->getEnergy()+this->electronMass));
                    double secondRatioOfEnergy = 1/firstRatioOfEnergy;
                    //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<" S3: "<<label3<<std::endl;
                    sumOfComponentA += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVector*spectralComponent3D);
                    sumOfComponentB += ((polarizationVector^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-((polarizationVector^(spectralComponent3D*this->getEnergy()-photonEmissionVector*this->getPhotonEnergy()*spectralComponent[0])*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))));
                }
            }
            std::complex<double> componentASpinCoefficient = (
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))                    )
            );
            Dimension3Vector<std::complex<double>> componentBSpinCoefficient = Dimension3Vector<std::complex<double>>(
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,1)
                )+
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))*
                    std::complex<double>(0,-1)
                ),
                (
                    (std::cos((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::cos((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,-(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )-
                (
                    (std::sin((axisOfEmissionPolarAngleOfElectronSpin/2.0)-((spinEmission-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfEmissionAzimuthalAngleOfElectronSpin/2))))*
                    (std::sin((axisOfIncidentPolarAngleOfElectronSpin/2.0)-((spinIncident-0.5)/2)*M_PI)*std::exp(std::complex<double>(0,(axisOfIncidentAzimuthalAngleOfElectronSpin/2))))
                )
            );
            sixTermsOfDifferentialEmissionIntensity0 += ((0.5+(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(std::conj(sumOfComponentA)*sumOfComponentA)).real();
            sixTermsOfDifferentialEmissionIntensity1 += ((0.5-(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(sumOfComponentB.conjugate()*sumOfComponentB)).real();
            sixTermsOfDifferentialEmissionIntensity2 += ((sumOfComponentB*combinedEmissionOrientationAxis)*(sumOfComponentB.conjugate()*combinedIncidentOrientationAxis)).real();
            sixTermsOfDifferentialEmissionIntensity3 += ((-(sumOfComponentB.conjugate()^sumOfComponentB)+(sumOfComponentB.conjugate()*sumOfComponentA)+(sumOfComponentB*std::conj(sumOfComponentA)))*(combinedEmissionOrientationAxis^combinedIncidentOrientationAxis)*(0.5)).real();
            sixTermsOfDifferentialEmissionIntensity4 += ((sumOfComponentB.conjugate()^sumOfComponentB)*(combinedEmissionOrientationAxis-combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)).real();
            sixTermsOfDifferentialEmissionIntensity5 += (((sumOfComponentB*(std::conj(sumOfComponentA)))-(sumOfComponentB.conjugate()*sumOfComponentA))*(combinedEmissionOrientationAxis+combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)).real();
        }
    }
    // theta+5/gamma
    //file.close();
    //std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    double fineStructureConstant = 1.0/137;
    return {sixTermsOfDifferentialEmissionIntensity0,sixTermsOfDifferentialEmissionIntensity1,sixTermsOfDifferentialEmissionIntensity2,sixTermsOfDifferentialEmissionIntensity3,sixTermsOfDifferentialEmissionIntensity4,sixTermsOfDifferentialEmissionIntensity5};
}

double RadiationWithSpinAndPolarzation::getSpinIncident() {
    return spinIncident;
}

double RadiationWithSpinAndPolarzation::getSpinEmission() {
    return spinEmission;
}

double RadiationWithSpinAndPolarzation::getAxisOfIncidentAzimuthalAngleOfElectronSpin() {
    return axisOfIncidentAzimuthalAngleOfElectronSpin;
}

double RadiationWithSpinAndPolarzation::getAxisOfIncidentPolarAngleOfElectronSpin() {
    return axisOfIncidentPolarAngleOfElectronSpin;
}

double RadiationWithSpinAndPolarzation::getAxisOfEmissionAzimuthalAngleOfElectronSpin() {
    return axisOfEmissionAzimuthalAngleOfElectronSpin;
}

double RadiationWithSpinAndPolarzation::getAxisOfEmissionPolarAngleOfElectronSpin() {
    return axisOfEmissionPolarAngleOfElectronSpin;
}

Dimension3Vector<double> RadiationWithSpinAndPolarzation::getIncidentOrientationAxis() {
    return incidentOrientationAxis;
}

Dimension3Vector<double> RadiationWithSpinAndPolarzation::getEmissionOrientationAxis() {
    return incidentOrientationAxis;
}

std::vector<double> RadiationWithSpinAndPolarzation::getStokesParameter() {
    return stokesParameter;
}

std::vector<double> RadiationWithSpinAndPolarzation::getStokesParameterNormalized() {
    return std::vector<double>({stokesParameter[0]/stokesParameter[0],stokesParameter[1]/stokesParameter[0],stokesParameter[2]/stokesParameter[0],stokesParameter[3]/stokesParameter[0]});
}

RadiationWithSpinAndPolarzation::~RadiationWithSpinAndPolarzation() {
}

