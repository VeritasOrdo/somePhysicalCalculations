#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<fstream>
#include<ctime>
#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../package/BasicMathFunctionDefinition/BasicMathFunctionDefinition.h"
//#include "../package/RadiationOfElectron/BasicRadiation/BasicRadiation.h"
#include "../package/RadiationOfElectron/RadiationWithSpinAndPolarzation/RadiationWithSpinAndPolarzation.h"
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286

/*int main(){
    double electronMass = 511000;
    double fieldParameter1 = 1e-300;
    double fieldParameter2 = 1.0;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 30000000000.0;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double omega = 1.55;
    double spinIncidient = -0.5;
    double spinEmission = 0.5;
    double polarizationAlpha = 0;
    double polarizationBeta = 0;
    double axisOfIncidentPolarAngleOfElectronSpin = 0;
    double axisOfIncidentAzimuthalAngleOfElectronSpin = 0;
    double axisOfEmissionPolarAngleOfElectronSpin = 0;
    double axisOfEmissionAzimuthalAngleOfElectronSpin = 0;
    double azimuthalAngleOfEmission = 0;
    //time
    std::clock_t start = std::clock();
    std::fstream file;
    file.open("downup+.txt",std::ios::out);
    for(double photonEnergyGeV=0.0;photonEnergyGeV<25;photonEnergyGeV+=1.0){
        if(photonEnergyGeV==0.0){
            file << photonEnergyGeV << "\t" << 0 << std::endl;
            continue;
        }
        double photonEnergy = photonEnergyGeV*1000000000;
        RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin);
        radiationWithSpinAndPolarzation.calculateDifferentialEmissionIntensity();
        file<<photonEnergyGeV<<"\t"<<2*M_PI*radiationWithSpinAndPolarzation.getDifferentialEmissionIntensity()<<std::endl;
    }
    file.close();
    std::clock_t end = std::clock();
    std::cout<<"Time: "<<(end-start)/(double)CLOCKS_PER_SEC<<"s"<<std::endl;
    return 0;
}*/

/*int main(){
    double electronMass = 511000;
    double fieldParameter1 = 20;
    double fieldParameter2 = 0.3;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 4*reducedMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double omega = 1.55;
    double spinIncidient = -0.5;
    double spinEmission = -0.5;
    double polarizationAlpha = 0;
    double polarizationBeta = 0;
    double axisOfIncidentPolarAngleOfElectronSpin = 0;
    double axisOfIncidentAzimuthalAngleOfElectronSpin = 0;
    double axisOfEmissionPolarAngleOfElectronSpin = 0;
    double axisOfEmissionAzimuthalAngleOfElectronSpin = 0;
    double azimuthalAngleOfEmission = 0;
    double photonEnergy = 0.0001*energyPrime;
    RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin);
    radiationWithSpinAndPolarzation.calculateDifferentialEmissionIntensity();
    std::cout << "differentialEmissionIntensity: " << 2*M_PI*radiationWithSpinAndPolarzation.getDifferentialEmissionIntensity() << std::endl;
    return 0;
}*/

/*int main(){
    double electronMass = 511000;
    double fieldParameter1 = 1e-300;
    double fieldParameter2 = 1;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 4*reducedMass;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double PI = 3.14159265358979323846;
    double omega = 1.55;
    double spinIncidient = 0.5;
    double spinEmission = 0.5;
    double polarizationAlpha = 0;
    double polarizationBeta = 0;
    double axisOfIncidentPolarAngleOfElectronSpin = 0;
    double axisOfIncidentAzimuthalAngleOfElectronSpin = 0;
    double axisOfEmissionPolarAngleOfElectronSpin = 0;
    double axisOfEmissionAzimuthalAngleOfElectronSpin = 0;
    double azimuthalAngleOfEmission = 0;
    double photonEnergy = 0.0001*energyPrime;
    RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin);
    radiationWithSpinAndPolarzation.calculateDifferentialEmissionIntensity();
    std::cout << "differentialEmissionIntensity: " << 2*M_PI*radiationWithSpinAndPolarzation.getDifferentialEmissionIntensity() << std::endl;
    return 0;
}*/

/*int main(){
    double electronMass = 511000;
    double fieldParameter1 = 1.6;
    double fieldParameter2 = 0.1;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 15330000;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double omega = 1.55;
    double spinIncidient = 0.5;
    double spinEmission = 0.5;
    double polarizationAlpha = 0;
    double polarizationBeta = 0;
    double axisOfIncidentPolarAngleOfElectronSpin = 0;
    double axisOfIncidentAzimuthalAngleOfElectronSpin = 0;
    double axisOfEmissionPolarAngleOfElectronSpin = 0;
    double axisOfEmissionAzimuthalAngleOfElectronSpin = 0;
    double azimuthalAngleOfEmission = 0;
    double photonEnergy = 1000;
    RadiationWithSpinAndPolarzation radiationWithSpinAndPolarzation(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,spinIncidient,spinEmission,polarizationAlpha,polarizationBeta,axisOfIncidentAzimuthalAngleOfElectronSpin,axisOfIncidentPolarAngleOfElectronSpin,axisOfEmissionAzimuthalAngleOfElectronSpin,axisOfEmissionPolarAngleOfElectronSpin,1,1,0,0.001);
    radiationWithSpinAndPolarzation.calculateDifferentialEmissionIntensity();
    //BasicRadiationOfElectronInCounterpropagatingLaser radiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,1,1);
    //radiationOfElectron.calculateDifferentialEmissionIntensity();
    //double radiation = radiationOfElectron.getDifferentialEmissionIntensity();
    double radiation = radiationWithSpinAndPolarzation.getDifferentialEmissionIntensity();
    std::cout << "differentialEmissionIntensity: " << 2*M_PI*radiation << std::endl;
    //std::cout << "cosPI/2: " << std::cos(M_PI/2) << std::endl;
    //std::cout << "cosPI/2_2: " << std::cos(PI/2) << std::endl;
    return 0;
}*/

/*int main(){
    double electronMass = 511000;
    double fieldParameter1 = 1.6;
    double fieldParameter2 = 0.1;
    double reducedMass = electronMass*std::sqrt(1+fieldParameter1*fieldParameter1+fieldParameter2*fieldParameter2);
    double energyPrime = 15330000;
    double momentumXPrime = 0;
    double momentumZPrime = std::sqrt(energyPrime*energyPrime-reducedMass*reducedMass-momentumXPrime*momentumXPrime);
    double omega = 1.55;
    double spinIncidient = 0.5;
    double spinEmission = 0.5;
    double polarizationAlpha = 0;
    double polarizationBeta = 0;
    double axisOfIncidentPolarAngleOfElectronSpin = 0;
    double axisOfIncidentAzimuthalAngleOfElectronSpin = 0;
    double axisOfEmissionPolarAngleOfElectronSpin = 0;
    double axisOfEmissionAzimuthalAngleOfElectronSpin = 0;
    double azimuthalAngleOfEmission = 0;
    double photonEnergy = 5517.4;
    BasicRadiationOfElectronInCounterpropagatingLaser radiationOfElectron(momentumZPrime,momentumXPrime,fieldParameter1,fieldParameter2,0,photonEnergy,azimuthalAngleOfEmission,1,1);
    radiationOfElectron.calculateEmissionIntensityAndPolarization();
    std::vector<std::map<std::pair<int,int>,std::pair<double,double>>> emissionIntensityAndPolarization = radiationOfElectron.getEmissionMapIntensityList();
    std::fstream file;
    std::cout<<"emissionIntensityAndPolarization.size(): "<<emissionIntensityAndPolarization.size()<<std::endl;
    file.open("emissionIntensityAndPolarization.txt",std::ios::out);
    for(int i=0;i<emissionIntensityAndPolarization.size();i++){
        std::cout<<"i: "<<i<<std::endl;
        for(auto it=emissionIntensityAndPolarization[i].begin();it!=emissionIntensityAndPolarization[i].end();it++){
            file << it->first.first << "\t" << it->first.second << "\t" << it->second.first << "\t" << it->second.second << std::endl;
        }
    }
    file.close();
    return 0;
}*/