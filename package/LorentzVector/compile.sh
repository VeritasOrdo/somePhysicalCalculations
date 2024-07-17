#!/bin/bash

# Compile the source code
g++ -c LorentzVector.cpp -std=c++17
ar -rsv libLorentzVector.a LorentzVector.o