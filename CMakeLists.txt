cmake_minimum_required(VERSION 3.12)

# Set the project name
project(Physics)

# Set the cxx standard to c++17
set(CMAKE_CXX_STANDARD 17)

# Add the target
add_executable(test ./calculation/test.cpp)

# Add the libraries
target_link_libraries(test PUBLIC 
    ${CMAKE_SOURCE_DIR}/package/RadiationOfElectron/BasicRadiation/libBasicRadiation.a 
    ${CMAKE_SOURCE_DIR}/package/ElectronInCounterpropagatingLaser/libElectronInCounterpropagatingLaser.a 
    ${CMAKE_SOURCE_DIR}/package/LorentzVector/libLorentzVector.a 
    ${CMAKE_SOURCE_DIR}/package/BasicMathFunctionDefinition/libBasicMathFunctionDefinition.a)
