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
    std::fstream parameterFile("CalculateStokesParameter.json5");
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
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double momentumXPrime = parameter["momentumXPrime"];
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double omega = parameter["omega"];
    double spinIncidient = parameter["incident"]["spin"];
    double spinEmission = parameter["emission"]["spin"];
    double polarizationAlpha = 0;
    double polarizationBeta = 0;
    double axisOfIncidentPolarAngleOfElectronSpin = parameter["incident"]["polarAngle"];
    double axisOfIncidentAzimuthalAngleOfElectronSpin = parameter["incident"]["azimuthalAngle"];
    double axisOfEmissionPolarAngleOfElectronSpin = parameter["emission"]["polarAngle"];
    double axisOfEmissionAzimuthalAngleOfElectronSpin = parameter["emission"]["azimuthalAngle"];
    double azimuthalAngleOfEmission = parameter["azimuthalAngleOfEmission"];
    double begin = parameter["plot"]["begin"];
    double end = parameter["plot"]["end"];
    int points = parameter["plot"]["points"];
    double step = (end-begin)/points;
    bool normalizeFlag = parameter["normalize"];
    std::string fileName = parameter["outputFile"];
    double emissionPolarAngleLimitMin = parameter["emissionPolarAngleLimit"]["min"];
    double emissionPolarAngleLimitMax = parameter["emissionPolarAngleLimit"]["max"];
    std::fstream file;
    file.open(fileName,std::ios::out);
    for(double photonEnergy=begin+step;photonEnergy<end;photonEnergy+=step){
        RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2,emissionPolarAngleLimitMin,emissionPolarAngleLimitMax);
        radiationWithSpinAndPolarzation.calculateStokesParameter();
        if(normalizeFlag){
            std::vector<double> stokesParameterNormalized = radiationWithSpinAndPolarzation.getStokesParameterNormalized();
            file<<photonEnergy<<"\t"<<stokesParameterNormalized[0]<<"\t"<<stokesParameterNormalized[1]<<"\t"<<stokesParameterNormalized[2]<<"\t"<<stokesParameterNormalized[3]<<std::endl;
        }
        else{
            std::vector<double> stokesParameter = radiationWithSpinAndPolarzation.getStokesParameter();
            file<<photonEnergy<<"\t"<<stokesParameter[0]<<"\t"<<stokesParameter[1]<<"\t"<<stokesParameter[2]<<"\t"<<stokesParameter[3]<<std::endl;
        }
    }
    file.close();
    //RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,1,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,rotationDirection1,rotationDirection2);
    //radiationWithSpinAndPolarzation.calculateStokesParameter();
    //std::vector<double> stokesParameter = radiationWithSpinAndPolarzation.getStokesParameterNormalized();
    //std::cout<<"Stokes Parameter: "<<stokesParameter[0]<<" "<<stokesParameter[1]<<" "<<stokesParameter[2]<<" "<<stokesParameter[3]<<std::endl;
    return 0;
}


