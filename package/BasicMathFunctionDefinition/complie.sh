#!/bin/bash

# Compile the source code
g++ -c BasicMathFunctionDefinition.cpp -L../package/gsl/lib -lgsl -lgslcblas -std=c++17
ar -rsv libBasicMathFunctionDefinition.a BasicMathFunctionDefinition.o