#include "NumericalSolutionOfElectronInCounterpropagatingLaser.h"
#include <iostream>

int func(double t, const double y[], double f[], void *params) {
    double *p = (double *)params;
    double omega = p[0];
    double fieldParameter1 = p[1];
    double fieldParameter2 = p[2];
    f[0] = y[1];
    f[1] = -fieldParameter1*omega*sin(omega*t-omega*f[5])+fieldParameter1*omega*f[6]*sin(omega*t-omega*f[5])-fieldParameter2*omega*sin(omega*t+omega*f[5])-fieldParameter2*omega*f[6]*sin(omega*t+omega*f[5]);
    f[2] = y[3];
    f[3] = fieldParameter1*omega*cos(omega*t-omega*f[5])-fieldParameter1*omega*f[6]*cos(omega*t-omega*f[5])+fieldParameter2*omega*cos(omega*t+omega*f[5])+fieldParameter2*omega*f[6]*cos(omega*t+omega*f[5]);
    f[4] = y[5];
    f[5] = fieldParameter1*omega*f[4]*cos(omega*t-omega*f[5])-fieldParameter1*omega*f[2]*sin(omega*t-omega*f[5])+fieldParameter2*omega*f[2]*sin(omega*t+omega*f[5])-fieldParameter2*omega*f[4]*cos(omega*t+omega*f[5]);
    return GSL_SUCCESS;
}

NumericalSolutionOfElectronInCounterpropagatingLaser::NumericalSolutionOfElectronInCounterpropagatingLaser(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTimeBegin,double properTimeEnd,double stepLength) {
    this->fieldParameter1 = fieldParameter1;
    this->fieldParameter2 = fieldParameter2;
    this->properTimeBegin = properTimeBegin;
    this->properTimeEnd = properTimeEnd;
    this->stepLength = stepLength;
    this->electronMomentumPrime = Dimension3Vector<double>(momentumXPrime,0,momentumZPrime);
    this->electronMomentum={};
    this->electronMomentum.push_back(this->electronMomentumPrime);
    this->electronCoordinatePrime = Dimension3Vector<double>(0,0,0);
    this->electronCoordinate={};
    this->electronCoordinate.push_back(this->electronCoordinatePrime);
    this->EnengyPrime = sqrt(this->electronMomentumPrime*this->electronMomentumPrime+this->electronMass*this->electronMass);
    this->Enengy={};
    this->Enengy.push_back(this->EnengyPrime);
    this->electronVelocityPrime = this->electronMomentumPrime/this->EnengyPrime;
    this->electronVelocity={};
    this->electronVelocity.push_back(this->electronVelocityPrime);
    this->properTime={};
    this->properTime.push_back(this->properTimeBegin);
    double y[6] = {this->electronCoordinatePrime.getX(),this->electronVelocityPrime.getX(),this->electronCoordinatePrime.getZ(),this->electronVelocityPrime.getZ(),this->electronCoordinatePrime.getY(),this->electronVelocityPrime.getY()};
    double params[3] = {this->omega,this->fieldParameter1,this->fieldParameter2};
    gsl_odeiv2_system sys = {func, NULL, 6, params};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    for (double t = this->properTimeBegin; t < this->properTimeEnd; t += this->stepLength) {
        int status = gsl_odeiv2_driver_apply(d, &t, t + this->stepLength, y);
        if (status != GSL_SUCCESS) {
            std::cout << "error, return value=" << status << std::endl;
            break;
        }
        this->electronCoordinate.push_back(Dimension3Vector<double>(y[0],0,y[2]));
        this->electronVelocity.push_back(Dimension3Vector<double>(y[1],0,y[3]));
        this->electronMomentum.push_back(Dimension3Vector<double>(this->electronVelocity.back()*this->Enengy.back()));
        this->Enengy.push_back(sqrt(this->electronMomentum.back()*this->electronMomentum.back()+this->electronMass*this->electronMass));
        this->properTime.push_back(t);
    }
    gsl_odeiv2_driver_free(d);
}

std::vector<Dimension3Vector<double>> NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronMomentum() {
    return this->electronMomentum;
}

Dimension3Vector<double> NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronMomentumPrime() {
    return this->electronMomentumPrime;
}

std::vector<Dimension3Vector<double>> NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronCoordinate() {
    return this->electronCoordinate;
}

Dimension3Vector<double> NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronCoordinatePrime() {
    return this->electronCoordinatePrime;
}

std::vector<Dimension3Vector<double>> NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronVelocity() {
    return this->electronVelocity;
}

Dimension3Vector<double> NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronVelocityPrime() {
    return this->electronVelocityPrime;
}

std::vector<double> NumericalSolutionOfElectronInCounterpropagatingLaser::getProperTime() {
    return this->properTime;
}

std::vector<double> NumericalSolutionOfElectronInCounterpropagatingLaser::getEnengy() {
    return this->Enengy;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getEnengyPrime() {
    return this->EnengyPrime;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getStepLength() {
    return this->stepLength;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getFieldParameter1() {
    return this->fieldParameter1;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getFieldParameter2() {
    return this->fieldParameter2;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getProperTimeBegin() {
    return this->properTimeBegin;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getProperTimeEnd() {
    return this->properTimeEnd;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getOmega() {
    return this->omega;
}

double NumericalSolutionOfElectronInCounterpropagatingLaser::getElectronMass() {
    return this->electronMass;
}

