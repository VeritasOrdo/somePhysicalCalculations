#!/bin/bash

# Compile the source code
g++ -c BasicRadiation.cpp -L../../LorentzVector -lLorentzVector -L../../ElectronInCounterpropagatingLaser -lElectronInCounterpropagatingLaser -L../../BasicMathFunctionDefinition -lBasicMathFunctionDefinition
ar -rsv libBasicRadiation.a BasicRadiation.o