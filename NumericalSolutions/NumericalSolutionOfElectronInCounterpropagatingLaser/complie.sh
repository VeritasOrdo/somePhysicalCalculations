#!/bin/bash

# Compile the source code
g++ -c NumericalSolutionOfElectronInCounterpropagatingLaser.cpp -L../../package/Dimension3Vector -lDimension3Vector -lgsl -lgslcblas
ar -rsv libNumericalSolutionOfElectronInCounterpropagatingLaser.a NumericalSolutionOfElectronInCounterpropagatingLaser.o