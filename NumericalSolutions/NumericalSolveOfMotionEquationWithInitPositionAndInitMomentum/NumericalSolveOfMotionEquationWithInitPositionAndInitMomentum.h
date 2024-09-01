#include "../NumericalSolveOfMotionEquationBase/NumericalSolveOfMotionEquationBase.h"

#pragma once

class NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum:public NumericalSolveOfMotionEquationBase
{
private:
    Dimension3Vector<double> initPosition;
    Dimension3Vector<double> initMomentum;
    double initEnergy;
public:
    NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum(double mass, Dimension3Vector<double> initPosition,Dimension3Vector<double> initMomentum,double timeBegin,double timeEnd,double stepLength,std::function<Dimension3Vector<double>(double t,Dimension3Vector<double> y,Dimension3Vector<double> dydt)> motionEquation);
    Dimension3Vector<double> callMotionEquation(double t,Dimension3Vector<double> y,Dimension3Vector<double> dydt);
    Dimension3Vector<double> getInitPosition();
    Dimension3Vector<double> getInitMomentum();
    double getMass();
    double getInitEnergy();
    std::vector<Dimension3Vector<double>> getPosition();
    std::vector<Dimension3Vector<double>> getVelocity();
    std::vector<Dimension3Vector<double>> getMomentum();
    std::vector<double> getEnergy();
    ~NumericalSolveOfMotionEquationWithInitPositionAndInitMomentum();
};