#!/bin/bash

# Compile the source code
g++ -c Dimension3Vector.cpp
ar -rsv libDimension3Vector.a Dimension3Vector.o