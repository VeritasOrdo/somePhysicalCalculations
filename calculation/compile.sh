#!/bin/bash

# Compile the source code
g++ -o test test.cpp -L../package/ElectronInCounterpropagatingLaser -lElectronInCounterpropagatingLaser -L../package/LorentzVector -lLorentzVector -L../package/BasicMathFunctionDefinition -lBasicMathFunctionDefinition -L../package/RadiationOfElectron/BasicRadiation -lBasicRadiation -L../package/BasicMathFunctionDefinition -lBasicMathFunctionDefinition -std=c++17
./test

