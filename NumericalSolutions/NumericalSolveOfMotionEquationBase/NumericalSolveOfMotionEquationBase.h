#include "../../package/Dimension3Vector/Dimension3Vector.h"
#include<cmath>
#include<vector>
#include<functional>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

#pragma once

class NumericalSolveOfMotionEquationBase
{
private:
    double timeBegin;
    double timeEnd;
    double stepLength;
protected:
    std::vector<double> time;
    std::function<Dimension3Vector<double>(double t,Dimension3Vector<double> y,Dimension3Vector<double> dydt)> motionEquation;
    std::vector<Dimension3Vector<double>> position;
    std::vector<Dimension3Vector<double>> velocity;
    std::vector<Dimension3Vector<double>> momentum; 
    std::vector<double> energy;
    double mass;
public:
    NumericalSolveOfMotionEquationBase(double mass, double timeBegin, double timeEnd, double stepLength, std::function<Dimension3Vector<double>(double t,Dimension3Vector<double> y,Dimension3Vector<double> dydt)> motionEquation);
    virtual std::vector<double> getTime();
    virtual double getMass();
    double getTimeBegin();
    double getTimeEnd();
    double getStepLength();
    ~NumericalSolveOfMotionEquationBase();
};

