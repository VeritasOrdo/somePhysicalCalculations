#include<fstream>
#include<iostream>
#include<chrono>
#include<nlohmann/json.hpp>
#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../package/BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
#include "../package/RadiationOfElectron/RadiationWithSpinAndPolarzation/RadiationWithSpinAndPolarzation.h"

using json = nlohmann::json;

int main(){
    std::fstream parameterFile("DrawRadiationPlots.json");
    if(!parameterFile.is_open()){
        std::cout<<"parameter file not found"<<std::endl;
        return 0;
    }
    json parameter= json::parse(parameterFile);
    double electronMass = 511000;
    double fieldParameter1 = parameter["fieldParameter1"];
    double fieldParameter2 = parameter["fieldParameter2"];
    double energyPrime = parameter["energyPrime"];
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double momentumXPrime = parameter["momentumXPrime"];
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
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
        RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin);
        radiationWithSpinAndPolarzation.calculateDifferentialEmissionIntensity();
        file<<photonEnergy<<"\t"<<2*M_PI*radiationWithSpinAndPolarzation.getDifferentialEmissionIntensity()<<std::endl;
    }
    file.close();
    //time end
    auto totalTimeEnd = std::chrono::system_clock::now();
    std::cout<<"Total time: "<<std::chrono::duration_cast<std::chrono::hours>(totalTimeEnd-totalTimeBegin).count()<<"h"<<std::chrono::duration_cast<std::chrono::minutes>(totalTimeEnd-totalTimeBegin).count()%60<<"m"<<std::chrono::duration_cast<std::chrono::seconds>(totalTimeEnd-totalTimeBegin).count()%60<<"s"<<std::endl;
    return 0;
}
