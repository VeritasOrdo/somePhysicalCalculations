#!/bin/bash

# Compile the source code
g++ -c BasicMathFunctionDefinition.cpp -L../gsl/include/g -lgsl -lgslcblas -std=c++17
ar -rsv libBasicMathFunctionDefinition.a BasicMathFunctionDefinition.o