#include "NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum.h"
#include <iostream>
#include <gsl/gsl_odeiv2.h> // Include the appropriate header file

int func(double t,const double y[],double f[],void *params) {
        //this->motionEquation(t,Dimension3Vector<double>(y[0],y[1],y[2]),Dimension3Vector<double>(y[3],y[4],y[5]));
        NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum *p = (NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum *)params;
        double energy = std::sqrt((p->getMass()*p->getMass())+((y[3]*y[3])+(y[4]*y[4])+(y[5]*y[5])));
        //std::cout<<"energy: "<<energy<<std::endl;
        f[0] = y[3]/energy;
        f[1] = y[4]/energy;
        f[2] = y[5]/energy;
        //std::cout<<"position: "<<Dimension3Vector<double>(y[0],y[1],y[2])<<std::endl;
        //std::cout<<"Momentum: "<<Dimension3Vector<double>(y[3],y[4],y[5])<<std::endl;
        Dimension3Vector<double> acceleration = p->callMotionEquation(t,Dimension3Vector<double>(y[0],y[1],y[2]),Dimension3Vector<double>(y[3],y[4],y[5]));
        //std::cout<<"acceleration: "<<acceleration<<std::endl;
        //std::cin.get();
        f[3] = acceleration.getX();
        f[4] = acceleration.getY();
        f[5] = acceleration.getZ();
        //std::cout<<"f[0]: "<<f[0]<<std::endl;
        return GSL_SUCCESS;
};

NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum(double mass, Dimension3Vector<double> initPosition, Dimension3Vector<double> InitMomentum, double timeBegin, double timeEnd, double stepLength, std::function<Dimension3Vector<double>(double t, Dimension3Vector<double> y, Dimension3Vector<double> dydt)> motionEquation) : NumericalSolveOfMotionEquationBase(mass, timeBegin, timeEnd, stepLength, motionEquation){
    this->initPosition = initPosition;
    this->initMomentum = InitMomentum;
    this->initEnergy = std::sqrt((this->getMass()*this->getMass())+(this->initMomentum*this->initMomentum));
    std::cout<<"initEnergy: "<<this->initEnergy<<std::endl;
    this->position.push_back(initPosition);
    this->momentum.push_back(initMomentum);
    this->velocity.push_back(this->initMomentum/this->initEnergy);
    std::cout<<"initVelocity: "<<this->velocity[0]<<std::endl;
    this->energy.push_back(this->initEnergy);
    double t = timeBegin;
    double y[6] = {initPosition.getX(), initPosition.getY(), initPosition.getZ(), initMomentum.getX(), initMomentum.getY(), initMomentum.getZ()};
    gsl_odeiv2_system sys = {func, NULL, 6, this};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, stepLength, 1e-6, 1e-6);
    for(int i=1;i<this->time.size();i++){
        //std::cout<<"========================="<<std::endl;
        if(i%1000==0){
            std::cout<<i<<"/"<<this->time.size()<<std::endl;
        }
        double ti = this->time[i];
        //std::cout<<"t: "<<t<<std::endl;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        //std::cout<<"t: "<<t<<std::endl;
        if (status != GSL_SUCCESS) {
            std::cout << "error, return value=" << status << std::endl;
            break;
        }
        this->position.push_back(Dimension3Vector<double>(y[0],y[1],y[2]));
        this->momentum.push_back(Dimension3Vector<double>(y[3],y[4],y[5]));
        this->energy.push_back(std::sqrt((this->getMass()*this->getMass())+(y[3]*y[3])+(y[4]*y[4])+(y[5]*y[5])));
        this->velocity.push_back(Dimension3Vector<double>(y[3],y[4],y[5])/this->energy[i]);
    }
    gsl_odeiv2_driver_free(d);
}

Dimension3Vector<double> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::callMotionEquation(double t, Dimension3Vector<double> y, Dimension3Vector<double> dydt) {
    return this->motionEquation(t, y, dydt);
}

Dimension3Vector<double> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getInitPosition() {
    return this->initPosition;
}

Dimension3Vector<double> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getInitMomentum() {
    return this->initMomentum;
}

double NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getMass() {
    return this->mass;
}

double NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getInitEnergy() {
    return this->initEnergy;
}

std::vector<Dimension3Vector<double>> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getPosition() {
    return this->position;
}

std::vector<Dimension3Vector<double>> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getVelocity() {
    return this->velocity;
}

std::vector<Dimension3Vector<double>> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getMomentum() {
    return this->momentum;
}

std::vector<double> NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::getEnergy() {
    return this->energy;
}

NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum::~NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum() {
    // Do nothing
}

