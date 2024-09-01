#!/bin/bash

# Compile the source code
g++ -o test test.cpp -L../NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum -lNumericalSolveOfMotionEquationWithInitPositionAndInitMomentum -L../NumericalSolveOfMotionEquationBase -lNumericalSolveOfMotionEquationBase -L../../package/Dimension3Vector -lDimension3Vector -lgsl -lgslcblas
./test

