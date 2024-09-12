#!/bin/bash

# Compile the source code
g++ -o CalculateStokesParameter.out CalculateStokesParameter.cpp -Wl,-Bstatic -L../package/RadiationOfElectron/RadiationWithSpinAndPolarzation -lRadiationWithSpinAndPolarzation -L../package/Dimension3Vector -lDimension3Vector  -L../package/RadiationOfElectron/BasicRadiation -lBasicRadiation -L../package/ElectronInCounterpropagatingLaser -lElectronInCounterpropagatingLaser -L../package/BasicMathFunctionDefinition -lBasicMathFunctionDefinition -L../package/LorentzVector -lLorentzVector -L../package/gsl/lib -lgsl -lgslcblas -fopenmp -Wl,-Bdynamic