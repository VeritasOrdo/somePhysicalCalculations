#!/bin/bash

# Compile the source code
g++ -o test test.cpp -L../NumericalSolutionOfElectronInCounterpropagatingLaser -lNumericalSolutionOfElectronInCounterpropagatingLaser -L../../package/Dimension3Vector -lDimension3Vector -lgsl -lgslcblas
./test

