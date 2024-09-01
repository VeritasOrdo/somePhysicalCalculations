#!/bin/bash

# Compile the source code
g++ -c NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum.cpp -L../NumericalSolveOfMotionEquationBase -lNumericalSolveBaseOfMotionEquation -L../../package/Dimension3Vector -lDimension3Vector -lgsl -lgslcblas
ar -rsv libNumericalSolveOfMotionEquationWithInitPositionAndInitMomentum.a NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum.o