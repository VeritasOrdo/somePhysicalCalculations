#!/bin/bash

# Compile the source code
g++ -o test test.cpp -L../package/RadiationOfElectron/BasicRadiation -lBasicRadiation -L../package/BasicMathFunctionDefinition -lBasicMathFunctionDefinition -L../package/ElectronInCounterpropagatingLaser -lElectronInCounterpropagatingLaser -L../package/LorentzVector -lLorentzVector -std=c++17
./test

