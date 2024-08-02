#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<fstream>
#include "../package/ElectronInCounterpropagatingLaser/ElectronInCounterpropagatingLaser.h"
#include "../package/RadiationOfElectron/BasicRadiation/BasicRadiation.h"
#include "../package/Dimension3Vector/Dimension3Vector.h"

int main(){
    Dimension3Vector<std::complex<double>> vector1(std::complex<double>(1,2),std::complex<double>(3,4),std::complex<double>(5,6));
    Dimension3Vector<std::complex<double>> vector2(std::complex<double>(7,8),std::complex<double>(9,10),std::complex<double>(11,12));
    Dimension3Vector<double> vector3(1,2,3);
    std::cout << Dimension3Vector<std::complex<double>>(vector3) << std::endl;
    std::cout << vector2*2.6 << std::endl;
    double norm = vector1.getLength();
    std::cout << norm << std::endl;
    return 0;
}

