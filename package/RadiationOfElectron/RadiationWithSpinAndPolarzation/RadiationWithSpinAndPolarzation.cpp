#include "RadiationWithSpinAndPolarzation.h"
#include <chrono>
#include <omp.h>
#include <fstream>  

RadiationWithSpinAndPolarzation::RadiationWithSpinAndPolarzation(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime,double photonEnergy,double emissionAzimuthalAngle,double spinIncident,double spinEmission,double polarizationAlpha,double polarizationBeta,double axisOfIncidentAzimuthalAngleOfElectronSpin,double axisOfIncidentPolarAngleOfElectronSpin,double axisOfEmissionAzimuthalAngleOfElectronSpin,double axisOfEmissionPolarAngleOfElectronSpin) : BasicRadiationOfElectronInCounterpropagatingLaser(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,properTime,photonEnergy,emissionAzimuthalAngle) {
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
}

void RadiationWithSpinAndPolarzation::calculateDifferentialEmissionIntensity() {
    std::cout<<"this is the function in RadiationWithSpinAndPolarzation"<<std::endl;
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    //labelLeftLimit = 10;
    //labelRightLimit = 30;
    label3Limit = 10;
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    std::fstream file;
    file.open("spectralComponent.txt", std::ios::out);
    double sumOfSpectralComponent = 0;
    double sumOfSpectralComponentImag = 0;
    double sumOfSpectralComponentTest = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponent) reduction(+:sumOfSpectralComponentImag) reduction(+:sumOfSpectralComponentTest) 
    for(int labelLeft = 0; labelLeft <=labelLeftLimit; labelLeft++) {
        if(labelLeft%100==0){
            std::cout<<"SL: "<<labelLeft<<std::endl;
        }
        for(int labelRight = -labelRightLimit; labelRight <=labelRightLimit; labelRight++) {
            if(labelRight%100==0){
                std::cout<<"SR: "<<labelRight<<std::endl;
            }
            std::vector<double> emissionPolarAngles = calculateEmissionPolarAngle(labelLeft,labelRight);
            std::complex<double> sumOfComponentA=0;
            std::complex<double> sumOfComponentAX=0;
            std::complex<double> sumOfComponentAY=0;
            Dimension3Vector<std::complex<double>> sumOfComponentB = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBX = Dimension3Vector<std::complex<double>>(0,0,0);
            Dimension3Vector<std::complex<double>> sumOfComponentBY = Dimension3Vector<std::complex<double>>(0,0,0);
            for(int label3 = -label3Limit; label3 <=label3Limit; label3++) {
                for(int labelOfEmissionPolarAngles = 0; labelOfEmissionPolarAngles < emissionPolarAngles.size(); labelOfEmissionPolarAngles++) {
                    std::vector<std::complex<double>> spectralComponent = SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngles[labelOfEmissionPolarAngles]);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1],spectralComponent[2],spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    //Dimension3Vector<std::complex<double>> spectralComponentPhotonEmissionVector = ((photonEmissionVector*spectralComponent[4])/this->electronMass)*this->getPhotonEnergy();
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVector = polarizationVectorBase1*std::cos(polarizationAlpha)+polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    polarizationVector = ((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized()+(((photonEmissionVector)^((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized())).normalized())*std::complex<double>(0,1))*(1.0/sqrt(2.0));
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
            std::complex<double> sumOfSpectralComponentAmplitude = sumOfComponentA*componentASpinCoefficient+sumOfComponentB*componentBSpinCoefficient*std::complex<double>(0,1);
            sumOfSpectralComponentTest += std::abs(sumOfSpectralComponentAmplitude)*std::abs(sumOfSpectralComponentAmplitude);
            sumOfSpectralComponent += (
                (0.5+(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(std::conj(sumOfComponentA)*sumOfComponentA)+
                (0.5-(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(sumOfComponentB.conjugate()*sumOfComponentB)+
                (sumOfComponentB*combinedEmissionOrientationAxis)*(sumOfComponentB.conjugate()*combinedIncidentOrientationAxis)-
                (-(sumOfComponentB.conjugate()^sumOfComponentB)+(sumOfComponentB.conjugate()*sumOfComponentA)+(sumOfComponentB*std::conj(sumOfComponentA)))*(combinedEmissionOrientationAxis^combinedIncidentOrientationAxis)*(0.5)+
                (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedEmissionOrientationAxis-combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)+
                ((sumOfComponentB*(std::conj(sumOfComponentA)))-(sumOfComponentB.conjugate()*sumOfComponentA))*(combinedEmissionOrientationAxis+combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)
            ).real();
            sumOfSpectralComponentImag += (
                (0.5+(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(std::conj(sumOfComponentA)*sumOfComponentA)+
                (0.5-(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(sumOfComponentB.conjugate()*sumOfComponentB)+
                (sumOfComponentB*combinedEmissionOrientationAxis)*(sumOfComponentB.conjugate()*combinedIncidentOrientationAxis)-
                (-(sumOfComponentB.conjugate()^sumOfComponentB)+(sumOfComponentB.conjugate()*sumOfComponentA)+(sumOfComponentB*std::conj(sumOfComponentA)))*(combinedEmissionOrientationAxis^combinedIncidentOrientationAxis)*(0.5)+
                (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedEmissionOrientationAxis-combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)+
                ((sumOfComponentB*(std::conj(sumOfComponentA)))-(sumOfComponentB.conjugate()*sumOfComponentA))*(combinedEmissionOrientationAxis+combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)
            ).imag();
            if(sumOfComponentA.real()!=0||sumOfComponentB!=Dimension3Vector<std::complex<double>>(0,0,0)){
                file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<std::endl;
                file<<(0.5+(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(std::conj(sumOfComponentA)*sumOfComponentA)<<std::endl;
                file<<(0.5-(0.5*(combinedEmissionOrientationAxis*combinedIncidentOrientationAxis)))*(sumOfComponentB.conjugate()*sumOfComponentB)<<std::endl;
                file<<(sumOfComponentB*combinedEmissionOrientationAxis)*(sumOfComponentB.conjugate()*combinedIncidentOrientationAxis)<<std::endl;
                file<<-(-(sumOfComponentB.conjugate()^sumOfComponentB)+(sumOfComponentB.conjugate()*sumOfComponentA)+(sumOfComponentB*std::conj(sumOfComponentA)))*(combinedEmissionOrientationAxis^combinedIncidentOrientationAxis)*(0.5)<<std::endl;
                file<<(sumOfComponentB.conjugate()^sumOfComponentB)*(combinedEmissionOrientationAxis-combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)<<std::endl;
                file<<((sumOfComponentB*(std::conj(sumOfComponentA)))-(sumOfComponentB.conjugate()*sumOfComponentA))*(combinedEmissionOrientationAxis+combinedIncidentOrientationAxis)*(0.5)*std::complex<double>(0,1)<<std::endl;
            }
        }
    }
    //file.close();
    std::cout<<"sumOfSpectralComponentReal: "<<sumOfSpectralComponent<<std::endl;
    std::cout<<"sumOfSpectralComponentImag: "<<sumOfSpectralComponentImag<<std::endl;
    std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    double fineStructureConstant = 1.0/137;
    this->setDifferentialEmissionIntensity(((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponent);
    std::cout<<"differentialEmissionIntensity: "<<this->getDifferentialEmissionIntensity()<<std::endl;
    //std::cout<<"differentialEmissionIntensityTest: "<<((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponentTest<<std::endl;
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

RadiationWithSpinAndPolarzation::~RadiationWithSpinAndPolarzation() {
}

