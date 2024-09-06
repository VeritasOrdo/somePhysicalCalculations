#!/bin/bash

# Compile the source code
g++ -c RadiationWithSpinAndPolarzation.cpp  -L../../BasicMathFunctionDefinition -lBasicMathFunctionDefinition -L../../ElectronInCounterpropagatingLaser -lElectronInCounterpropagatingLaser -L../../LorentzVector -lLorentzVector -L../package/gsl/lib -lgsl -lgslcblas -fopenmp -std=c++17
ar -rsv libRadiationWithSpinAndPolarzation.a RadiationWithSpinAndPolarzation.o