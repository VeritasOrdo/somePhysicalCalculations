#include "RadiationWithSpinAndPolarzation.h"
#include <chrono>

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
}

void RadiationWithSpinAndPolarzation::calculateDifferentialEmissionIntensity() {
    std::vector<int> labelLimits = calculateLabelLimits();
    int labelLeftLimit = labelLimits[0]+labelLimits[2];
    int labelRightLimit = labelLimits[1]+labelLimits[2];
    int label3Limit = labelLimits[2];
    std::cout << "labelLeft: " << labelLeftLimit << std::endl;
    std::cout << "labelRight: " << labelRightLimit << std::endl;
    std::cout << "label3: " << label3Limit << std::endl;
    double sumOfSpectralComponent = 0;
    long double time0 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    #pragma omp parallel for schedule(dynamic) reduction(+:sumOfSpectralComponent)
    for(int labelLeft = 0; labelLeft <=labelLeftLimit; labelLeft++) {
        for(int labelRight = -labelRightLimit; labelRight <=labelRightLimit; labelRight++) {
            std::vector<double> emissionPolarAngles = calculateEmissionPolarAngle(labelLeft,labelRight);
            std::complex<double> sumOfComponentA=0;
            Dimension3Vector<std::complex<double>> sumOfComponentB = Dimension3Vector<std::complex<double>>(0,0,0);
            for(int label3 = -label3Limit; label3 <=label3Limit; label3++) {
                for(int labelOfEmissionPolarAngles = 0; labelOfEmissionPolarAngles < emissionPolarAngles.size(); labelOfEmissionPolarAngles++) {
                    std::vector<std::complex<double>> spectralComponent = SpectralComponent(labelLeft,labelRight,label3,emissionPolarAngles[labelOfEmissionPolarAngles]);
                    Dimension3Vector<std::complex<double>> spectralComponent3D = Dimension3Vector<std::complex<double>>(spectralComponent[1],spectralComponent[2],spectralComponent[3]);
                    Dimension3Vector<std::complex<double>> photonEmissionVector = Dimension3Vector<std::complex<double>>(std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::sin(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> spectralComponentPhotonEmissionVector = photonEmissionVector*(spectralComponent[4]/this->electronMass);
                    Dimension3Vector<std::complex<double>> polarizationVectorBase1 = Dimension3Vector<std::complex<double>>(std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::cos(this->getEmissionAzimuthalAngle()),std::cos(emissionPolarAngles[labelOfEmissionPolarAngles])*std::sin(this->getEmissionAzimuthalAngle()),-std::sin(emissionPolarAngles[labelOfEmissionPolarAngles]));
                    Dimension3Vector<std::complex<double>> polarizationVectorBase2 = Dimension3Vector<std::complex<double>>(-std::sin(this->getEmissionAzimuthalAngle()),std::cos(this->getEmissionAzimuthalAngle()),0);
                    Dimension3Vector<std::complex<double>> polarizationVector = polarizationVectorBase1*std::cos(polarizationAlpha)+polarizationVectorBase2*std::sin(polarizationAlpha)*std::exp(std::complex<double>(0,polarizationBeta));
                    Dimension3Vector<std::complex<double>> polarizationVectorConjugate = polarizationVector.conjugate();
                    double firstRatioOfEnergy = std::sqrt((this->getResidualEnergy()+this->electronMass)/(this->getEnergy()+this->electronMass));
                    double secondRatioOfEnergy = 1/firstRatioOfEnergy;
                    sumOfComponentA += ((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())))*(firstRatioOfEnergy+secondRatioOfEnergy)*(polarizationVectorConjugate*spectralComponent3D);
                    sumOfComponentB += ((polarizationVectorConjugate^spectralComponent3D)*firstRatioOfEnergy-(polarizationVectorConjugate^(spectralComponent3D-photonEmissionVector)*secondRatioOfEnergy))*((this->getEnergy())/(2*std::sqrt(this->getEnergy()*this->getResidualEnergy())));
                }
            }
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
    long double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Time: " << time1-time0 << std::endl;
    double fineStructureConstant = 1.0/137;
    this->setDifferentialEmissionIntensity(((fineStructureConstant*electronMass*electronMass*this->getPhotonEnergy()*this->getResidualEnergy())/(2.0*M_PI*this->getEnergy()*this->getEnergy()*this->getEnergy()))*sumOfSpectralComponent);
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

