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
    this->incidentOrientationAxis = Dimension3Vector<double>(cos(axisOfIncidentAzimuthalAngleOfElectronSpin)*sin(axisOfIncidentPolarAngleOfElectronSpin),sin(axisOfIncidentAzimuthalAngleOfElectronSpin)*sin(axisOfIncidentPolarAngleOfElectronSpin),cos(axisOfIncidentPolarAngleOfElectronSpin));
    this->emissionOrientationAxis = Dimension3Vector<double>(cos(axisOfEmissionAzimuthalAngleOfElectronSpin)*sin(axisOfEmissionPolarAngleOfElectronSpin),sin(axisOfEmissionAzimuthalAngleOfElectronSpin)*sin(axisOfEmissionPolarAngleOfElectronSpin),cos(axisOfEmissionPolarAngleOfElectronSpin));
    this->combinedIncidentOrientationAxis = incidentOrientationAxis*this->spinIncident*2.0;
    this->combinedEmissionOrientationAxis = emissionOrientationAxis*this->spinEmission*2.0;
    this->stokesParameterNormalized = {};
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
    //std::fstream file;
    //file.open("spectralComponent.txt", std::ios::out);
    double sumOfSpectralComponent = 0;
    double sumOfSpectralComponentImag = 0;
    double sumOfSpectralComponentTest = 0;
    double sumOfSpectralComponentX = 0;
    double sumOfSpectralComponentY = 0;
    double stokesParameterI = 0;
    double stokesParameterQ = 0;
    double stokesParameterU = 0;
    double stokesParameterV = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponent) reduction(+:sumOfSpectralComponentTest) reduction(+:sumOfSpectralComponentX) reduction(+:sumOfSpectralComponentY) reduction(+:stokesParameterI) reduction(+:stokesParameterQ) reduction(+:stokesParameterU) reduction(+:stokesParameterV)
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
                    Dimension3Vector<std::complex<double>> spectralComponentPhotonEmissionVector = ((photonEmissionVector*spectralComponent[4])/this->electronMass)*this->getPhotonEnergy();
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVector = polarizationVectorBase1*std::cos(polarizationAlpha)+polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    polarizationVector = ((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized()+(((photonEmissionVector)^((Dimension3Vector<std::complex<double>>(0,1,0)^photonEmissionVector).normalized())).normalized())*std::complex<double>(0,1))*(1.0/sqrt(2.0));
                    Dimension3Vector<std::complex<double>> polarizationVectorX = polarizationVectorBase1*std::cos(polarizationAlpha);
                    Dimension3Vector<std::complex<double>> polarizationVectorY = polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    polarizationVectorX = polarizationVectorBase1*(polarizationVector*polarizationVectorBase1);
                    polarizationVectorY = polarizationVectorBase2*(polarizationVector*polarizationVectorBase2);
                    Dimension3Vector<std::complex<double>> polarizationVectorConjugate = polarizationVector.conjugate();
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy()+this->electronMass)/(this->getEnergy()+this->electronMass));
                    double secondRatioOfEnergy = 1/firstRatioOfEnergy;
                    //file<<"SL: "<<labelLeft<<" SR: "<<labelRight<<" S3: "<<label3<<std::endl;
                    sumOfComponentA += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorConjugate*spectralComponent3D);
                    sumOfComponentAX += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorX.conjugate()*spectralComponent3D);
                    sumOfComponentAY += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorY.conjugate()*spectralComponent3D);
                    //sumOfComponentB += (polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy()/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-(polarizationVectorConjugate^(spectralComponent3D*this->getEnergy()-photonEmissionVector*spectralComponent[0]))*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    //sumOfComponentBX += (polarizationVectorX.conjugate()^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy()/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-(polarizationVectorX.conjugate()^(spectralComponent3D*this->getEnergy()-photonEmissionVector*spectralComponent[0]))*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    //sumOfComponentBY += (polarizationVectorY.conjugate()^spectralComponent3D)*firstRatioOfEnergy*(this->getEnergy()/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))-(polarizationVectorY.conjugate()^(spectralComponent3D*this->getEnergy()-photonEmissionVector*spectralComponent[0]))*secondRatioOfEnergy*(1.0/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    /*if(spectralComponent3D!=Dimension3Vector<std::complex<double>>((0,0),(0,0),(0,0))){
                        file <<polarizationVectorConjugate<<'\t'<<spectralComponent3D<<std::endl;
                    }*/
                    sumOfComponentB += ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorConjugate^(spectralComponent3D-spectralComponentPhotonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    /*file<<"=============="<<std::endl;
                    if(spectralComponentPhotonEmissionVector!=Dimension3Vector<std::complex<double>>((0,0),(0,0),(0,0))){
                        file << polarizationVectorConjugate << std::endl;
                        file << spectralComponent3D << '\t' << spectralComponentPhotonEmissionVector << std::endl;
                        file << (polarizationVectorConjugate^spectralComponent3D) <<'\t' << (polarizationVectorConjugate^(spectralComponent3D-spectralComponentPhotonEmissionVector)) << std::endl;
                        file << ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorConjugate^(spectralComponent3D-spectralComponentPhotonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy()))) << std::endl;
                    }*/
                    sumOfComponentBX += ((polarizationVectorX.conjugate()^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorX.conjugate()^(spectralComponent3D-spectralComponentPhotonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                    sumOfComponentBY += ((polarizationVectorY.conjugate()^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorY.conjugate()^(spectralComponent3D-spectralComponentPhotonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                }
            }
            //file << "SL: " << labelLeft << " SR: " << labelRight << std::endl;
            //file << "sumOfComponentA: " << sumOfComponentA << std::endl;
            //file << "sumOfComponentB: " << sumOfComponentB << std::endl;
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
            //file<<sumOfComponentA*componentASpinCoefficient+sumOfComponentB*componentBSpinCoefficient<<std::endl;
            //file<<sumOfComponentA*componentASpinCoefficient-sumOfComponentB*componentBSpinCoefficient<<std::endl;
            std::complex<double> sumOfSpectralComponentAmplitudeX = sumOfComponentAX*componentASpinCoefficient+sumOfComponentBX*componentBSpinCoefficient*std::complex<double>(0,1);
            std::complex<double> sumOfSpectralComponentAmplitudeY = sumOfComponentAY*componentASpinCoefficient+sumOfComponentBY*componentBSpinCoefficient*std::complex<double>(0,1);
            stokesParameterI += std::abs(sumOfSpectralComponentAmplitudeX)*std::abs(sumOfSpectralComponentAmplitudeX)+std::abs(sumOfSpectralComponentAmplitudeY)*std::abs(sumOfSpectralComponentAmplitudeY);
            stokesParameterQ += std::abs(sumOfSpectralComponentAmplitudeX)*std::abs(sumOfSpectralComponentAmplitudeX)-std::abs(sumOfSpectralComponentAmplitudeY)*std::abs(sumOfSpectralComponentAmplitudeY);
            stokesParameterU += 2.0*std::real(sumOfSpectralComponentAmplitudeX*std::conj(sumOfSpectralComponentAmplitudeY));
            stokesParameterV += -2.0*std::imag(sumOfSpectralComponentAmplitudeX*std::conj(sumOfSpectralComponentAmplitudeY));
            sumOfSpectralComponentX += std::abs(sumOfSpectralComponentAmplitudeX)*std::abs(sumOfSpectralComponentAmplitudeX);
            sumOfSpectralComponentY += std::abs(sumOfSpectralComponentAmplitudeY)*std::abs(sumOfSpectralComponentAmplitudeY);
            sumOfSpectralComponentTest += std::abs(sumOfSpectralComponentAmplitude)*std::abs(sumOfSpectralComponentAmplitude);
            /*file<<(0.5+0.5*(combinedIncidentOrientationAxis*combinedEmissionOrientationAxis))*(std::conj(sumOfComponentA)*sumOfComponentA)<<'\t'<<
            (0.5-0.5*(combinedIncidentOrientationAxis*combinedEmissionOrientationAxis))*(sumOfComponentB.conjugate()*sumOfComponentB)<<'\t'<<
            (sumOfComponentB*combinedIncidentOrientationAxis)*(sumOfComponentB.conjugate()*combinedEmissionOrientationAxis)<<'\t'<<
            (-sumOfComponentB.conjugate()^sumOfComponentB+sumOfComponentB.conjugate()*sumOfComponentA+sumOfComponentB*std::conj(sumOfComponentA))*(combinedIncidentOrientationAxis^combinedEmissionOrientationAxis)*(-0.5)<<'\t'<<
            (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedIncidentOrientationAxis-combinedEmissionOrientationAxis)*(0.5)*std::complex<double>(0,1)<<'\t'<<
            (sumOfComponentB*(std::conj(sumOfComponentA))-sumOfComponentB.conjugate()*sumOfComponentA)*(combinedIncidentOrientationAxis+combinedEmissionOrientationAxis)*(0.5)*std::complex<double>(0,1)<<std::endl; 
            file<<(
                (0.5+0.5*(combinedIncidentOrientationAxis*combinedEmissionOrientationAxis))*(std::conj(sumOfComponentA)*sumOfComponentA)+
                (0.5-0.5*(combinedIncidentOrientationAxis*combinedEmissionOrientationAxis))*(sumOfComponentB.conjugate()*sumOfComponentB)+
                (sumOfComponentB*combinedIncidentOrientationAxis)*(sumOfComponentB.conjugate()*combinedEmissionOrientationAxis)+
                (-sumOfComponentB.conjugate()^sumOfComponentB+sumOfComponentB.conjugate()*sumOfComponentA+sumOfComponentB*std::conj(sumOfComponentA))*(combinedIncidentOrientationAxis^combinedEmissionOrientationAxis)*(-0.5)+
                (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedIncidentOrientationAxis-combinedEmissionOrientationAxis)*(0.5)*std::complex<double>(0,1)+
                (sumOfComponentB*(std::conj(sumOfComponentA))-sumOfComponentB.conjugate()*sumOfComponentA)*(combinedIncidentOrientationAxis+combinedEmissionOrientationAxis)*(0.5)*std::complex<double>(0,1)
            )<<std::endl;*/
            sumOfSpectralComponent += (
                (0.5+0.5*(combinedIncidentOrientationAxis*combinedEmissionOrientationAxis))*(std::conj(sumOfComponentA)*sumOfComponentA)+
                (0.5-0.5*(combinedIncidentOrientationAxis*combinedEmissionOrientationAxis))*(sumOfComponentB.conjugate()*sumOfComponentB)+
                (sumOfComponentB*combinedIncidentOrientationAxis)*(sumOfComponentB.conjugate()*combinedEmissionOrientationAxis)+
                (-sumOfComponentB.conjugate()^sumOfComponentB+sumOfComponentB.conjugate()*sumOfComponentA+sumOfComponentB*std::conj(sumOfComponentA))*(combinedIncidentOrientationAxis^combinedEmissionOrientationAxis)*(-0.5)+
                (sumOfComponentB.conjugate()^sumOfComponentB)*(combinedIncidentOrientationAxis-combinedEmissionOrientationAxis)*(0.5)*std::complex<double>(0,1)+
                (sumOfComponentB*(std::conj(sumOfComponentA))-sumOfComponentB.conjugate()*sumOfComponentA)*(combinedIncidentOrientationAxis+combinedEmissionOrientationAxis)*(0.5)*std::complex<double>(0,1)
            ).real();
        }
    }
    //file.close();
    std::cout<<"stokesParameterI: "<<stokesParameterI<<std::endl;
    std::cout<<"stokesParameterQ: "<<stokesParameterQ<<std::endl;
    std::cout<<"stokesParameterU: "<<stokesParameterU<<std::endl;
    std::cout<<"stokesParameterV: "<<stokesParameterV<<std::endl;
    this->stokesParameterNormalized.push_back(stokesParameterI/stokesParameterI);
    this->stokesParameterNormalized.push_back(stokesParameterQ/stokesParameterI);
    this->stokesParameterNormalized.push_back(stokesParameterU/stokesParameterI);
    this->stokesParameterNormalized.push_back(stokesParameterV/stokesParameterI);
    std::cout<<"stokesParameterNormalized: "<<"("<<this->stokesParameterNormalized[0]<<","<<this->stokesParameterNormalized[1]<<","<<this->stokesParameterNormalized[2]<<","<<this->stokesParameterNormalized[3]<<")"<<std::endl;
    std::cout<<"sumOfSpectralComponent: "<<sumOfSpectralComponent<<std::endl;
    std::cout<<"sumOfSpectralComponentImag: "<<sumOfSpectralComponentImag<<std::endl;
    std::cout<<"sumOfSpectralComponentTest: "<<sumOfSpectralComponentTest<<std::endl;
    std::cout<<"sumOfSpectralComponentX: "<<sumOfSpectralComponentX<<std::endl;
    std::cout<<"sumOfSpectralComponentY: "<<sumOfSpectralComponentY<<std::endl;
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << (time1-time0)/1000000000 << " seconds" << std::endl;
    double fineStructureConstant = 1.0/137;
    this->setDifferentialEmissionIntensity(((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponent);
    //std::cout<<"differentialEmissionIntensity: "<<this->getDifferentialEmissionIntensity()<<std::endl;
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

std::vector<double> RadiationWithSpinAndPolarzation::getStokesParameterNormalized() {
    return stokesParameterNormalized;
}

RadiationWithSpinAndPolarzation::~RadiationWithSpinAndPolarzation() {
}

