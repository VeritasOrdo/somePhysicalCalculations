#include<fstream>
#include<iostream>
#include<chrono>
#include<nlohmann/json.hpp>
#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../package/BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
#include "../package/RadiationOfElectron/BasicRadiation/BasicRadiation.h"
#include "../package/RadiationOfElectron/RadiationWithSpinAndPolarzation/RadiationWithSpinAndPolarzation.h"

using json = nlohmann::json;

int main(){
    std::fstream parameterFile("DrawRadiationPlots.json5");
    if(!parameterFile.is_open()){
        std::cout<<"parameter file not found"<<std::endl;
        return 0;
    }
    json parameter= json::parse(parameterFile);
    double electronMass = 511000;
    double fieldParameter1 = parameter["fieldParameter1"];
    double fieldParameter2 = parameter["fieldParameter2"];
    double rotationDirection1 = parameter["rotationDirection1"];
    double rotationDirection2 = parameter["rotationDirection2"];
    double energyPrime = parameter["energyPrime"];
    double reducedMass = electronMass*std::sqrt(1+(fieldParameter1*fieldParameter1)+(fieldParameter2*fieldParameter2));
    std::cout<<"reducedMass: "<<reducedMass<<std::endl;
    double momentumXPrime = parameter["momentumXPrime"];
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    std::cout<<"momentumZPrime: "<<momentumZPrime<<std::endl;
    std::cin.get();
    double omega = parameter["omega"];
    double spinIncidient = parameter["incident"]["spin"];
    double spinEmission = parameter["emission"]["spin"];
    double polarizationAlpha = parameter["polarization"]["alpha"];
    double polarizationBeta = parameter["polarization"]["beta"];
    double axisOfIncidentPolarAngleOfElectronSpin = parameter["incident"]["polarAngle"];
    double axisOfIncidentAzimuthalAngleOfElectronSpin = parameter["incident"]["azimuthalAngle"];
    double axisOfEmissionPolarAngleOfElectronSpin = parameter["emission"]["polarAngle"];
    double axisOfEmissionAzimuthalAngleOfElectronSpin = parameter["emission"]["azimuthalAngle"];
    double azimuthalAngleOfEmission = parameter["azimuthalAngleOfEmission"];
    double begin = parameter["plot"]["begin"];
    double end = parameter["plot"]["end"];
    double step = (end-begin)/(double(parameter["plot"]["points"]));
    bool useIncidentSpin = parameter["useIncidentSpin"];
    bool useEmissionSpin = parameter["useEmissionSpin"];
    bool usePolarization = parameter["usePolarization"];
    std::string fileName = parameter["outputFile"];
    //time start
    auto totalTimeBegin = std::chrono::system_clock::now();
    std::fstream file;
    file.open(fileName,std::ios::out);
    for(double photonEnergy=begin;photonEnergy<end;photonEnergy+=step){
        if(photonEnergy==0.0){
            file << photonEnergy << "\t" << 0 << std::endl;
            continue;
        }
        if((!useIncidentSpin)&&(!useEmissionSpin)&&(!usePolarization)){
            BasicRadiationOfElectronInCounterpropagatingLaser radiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,rotationDirection1,rotationDirection2);
            radiationOfElectron.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<2*M_PI*radiationOfElectron.getDifferentialEmissionIntensity()<<std::endl;
            continue;
        }
        if(useIncidentSpin&&useEmissionSpin&&(!usePolarization)){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation1(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,0,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation2(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,M_PI/2,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            //radiationWithSpinAndPolarzation1.calculateDifferentialEmissionIntensity();
            //radiationWithSpinAndPolarzation2.calculateDifferentialEmissionIntensity();
            std::vector<double> six1 = radiationWithSpinAndPolarzation1.calculateSixTermsOfDifferentialEmissionIntensity();
            std::vector<double> six2 = radiationWithSpinAndPolarzation2.calculateSixTermsOfDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<2*M_PI*(six1[0]+six2[0])<<'\t'<<2*M_PI*(six1[1]+six2[1])<<'\t'<<2*M_PI*(six1[2]+six2[2])<<'\t'<<2*M_PI*(six1[3]+six2[3])<<'\t'<<2*M_PI*(six1[4]+six2[4])<<'\t'<<2*M_PI*(six1[5]+six2[5])<<std::endl;
        }
        if((!useIncidentSpin)&&useEmissionSpin&&usePolarization){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation1(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,0.5,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation2(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,-0.5,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            radiationWithSpinAndPolarzation1.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation2.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<M_PI*(radiationWithSpinAndPolarzation1.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation2.getDifferentialEmissionIntensity())<<std::endl;
            continue;
        }
        if(useIncidentSpin&&(!useEmissionSpin)&&usePolarization){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation1(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,0.5,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation2(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,-0.5,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            radiationWithSpinAndPolarzation1.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation2.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<2*M_PI*(radiationWithSpinAndPolarzation1.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation2.getDifferentialEmissionIntensity())<<std::endl;
            continue;
        }
        if((!useIncidentSpin)&&(!useEmissionSpin)&&usePolarization){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation1(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,0.5,0.5,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation2(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,0.5,-0.5,polarizationAlpha+M_PI/2,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation3(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,-0.5,0.5,polarizationAlpha,polarizationBeta+M_PI/2,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation4(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,-0.5,-0.5,polarizationAlpha+M_PI/2,polarizationBeta+M_PI/2,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            radiationWithSpinAndPolarzation1.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation2.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation3.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation4.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<M_PI*(radiationWithSpinAndPolarzation1.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation2.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation3.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation4.getDifferentialEmissionIntensity())<<std::endl;
            continue;
        }
        if(useIncidentSpin&&(!useEmissionSpin)&&(!usePolarization)){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation1(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,0.5,0,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation2(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,-0.5,0,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation3(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,0.5,M_PI/2,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation4(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,-0.5,M_PI/2,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            radiationWithSpinAndPolarzation1.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation2.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation3.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation4.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<2*M_PI*(radiationWithSpinAndPolarzation1.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation2.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation3.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation4.getDifferentialEmissionIntensity())<<std::endl;
            continue;
        }
        if((!useIncidentSpin)&&useEmissionSpin&&(!usePolarization)){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation1(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,0.5,spinEmission,0,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation2(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,-0.5,spinEmission,0,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation3(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,0.5,spinEmission,M_PI/2,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation4(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,-0.5,spinEmission,M_PI/2,0,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            radiationWithSpinAndPolarzation1.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation2.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation3.calculateDifferentialEmissionIntensity();
            radiationWithSpinAndPolarzation4.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<M_PI*(radiationWithSpinAndPolarzation1.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation2.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation3.getDifferentialEmissionIntensity()+radiationWithSpinAndPolarzation4.getDifferentialEmissionIntensity())<<std::endl;
            continue;
        }
        if(useIncidentSpin&&useEmissionSpin&&usePolarization){
            RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
            radiationWithSpinAndPolarzation.calculateDifferentialEmissionIntensity();
            file<<photonEnergy<<"\t"<<2*M_PI*radiationWithSpinAndPolarzation.getDifferentialEmissionIntensity()<<std::endl;
            continue;
        }
    }
    file.close();
    //time end
    auto totalTimeEnd = std::chrono::system_clock::now();
    std::cout<<"Total time: "<<std::chrono::duration_cast<std::chrono::hours>(totalTimeEnd-totalTimeBegin).count()<<"h"<<std::chrono::duration_cast<std::chrono::minutes>(totalTimeEnd-totalTimeBegin).count()%60<<"m"<<std::chrono::duration_cast<std::chrono::seconds>(totalTimeEnd-totalTimeBegin).count()%60<<"s"<<std::endl;
    return 0;
}