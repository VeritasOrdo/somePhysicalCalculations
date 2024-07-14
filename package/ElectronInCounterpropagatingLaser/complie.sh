#!/bin/bash

# Compile the source code
g++ -c ElectronInCounterpropagatingLaser.cpp -L../LorentzVector -lLorentzVector
ar -rsv libElectronInCounterpropagatingLaser.a ElectronInCounterpropagatingLaser.o