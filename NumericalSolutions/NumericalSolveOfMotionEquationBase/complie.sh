#!/bin/bash

# Compile the source code
g++ -c NumericalSolveOfMotionEquationBase.cpp -L../../package/Dimension3Vector -lDimension3Vector -lgsl -lgslcblas
ar -rsv libNumericalSolveOfMotionEquationBase.a NumericalSolveOfMotionEquationBase.o