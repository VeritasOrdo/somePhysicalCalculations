#!/bin/bash

# Compile the source code
g++ -c BasicMathFunctionDefinition.cpp -std=c++17
ar -rsv libBasicMathFunctionDefinition.a BasicMathFunctionDefinition.o