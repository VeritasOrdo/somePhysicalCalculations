cmake_minimum_required(VERSION 3.12)

# Set the project name
project(Physics)

# Set the cxx standard to c++17
set(CMAKE_CXX_STANDARD 17)

# Add the target
add_executable(test ./calculation/test.cpp)

# Add the libraries
target_link_libraries(test PUBLIC 
    ${CMAKE_SOURCE_DIR}/package/RadiationOfElectron/RadiationWithSpinAndPolarzation/libRadiationWithSpinAndPolarzation.a
    ${CMAKE_SOURCE_DIR}/package/RadiationOfElectron/BasicRadiation/libBasicRadiation.a 
    ${CMAKE_SOURCE_DIR}/package/ElectronInCounterpropagatingLaser/libElectronInCounterpropagatingLaser.a 
    ${CMAKE_SOURCE_DIR}/package/LorentzVector/libLorentzVector.a 
    ${CMAKE_SOURCE_DIR}/package/Dimension3Vector/libDimension3Vector.a
    ${CMAKE_SOURCE_DIR}/package/BasicMathFunctionDefinition/libBasicMathFunctionDefinition.a
    ${CMAKE_SOURCE_DIR}/../gsl/lib/libgsl.a
    ${CMAKE_SOURCE_DIR}/../gsl/lib/libgslcblas.a)

# find the openmp
find_package(OpenMP REQUIRED)

# link the openmp
if(OpenMP_CXX_FOUND)
    target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)    
endif(OpenMP_CXX_FOUND)
