#include "NumericalSolveOfMotionEquationBase.h"

NumericalSolveOfMotionEquationBase::NumericalSolveOfMotionEquationBase(double mass, double timeBegin, double timeEnd, double stepLength, std::function<Dimension3Vector<double>(double t, Dimension3Vector<double> y, Dimension3Vector<double> dydt)> motionEquation)
{
    this->mass = mass;
    this->timeBegin = timeBegin;
    this->timeEnd = timeEnd;
    this->stepLength = stepLength;
    this->motionEquation = motionEquation;
    this->position = {};
    this->velocity = {};
    this->momentum = {};
    this->energy = {};
    this->time = {};
    this->time.push_back(timeBegin);
    for(double t = timeBegin; t < timeEnd; t += stepLength)
    {
        this->time.push_back(t + stepLength);
    }
}

double NumericalSolveOfMotionEquationBase::getMass()
{
    return this->mass;
}

std::vector<double> NumericalSolveOfMotionEquationBase::getTime()
{
    return this->time;
}

double NumericalSolveOfMotionEquationBase::getTimeBegin()
{
    return this->timeBegin;
}

double NumericalSolveOfMotionEquationBase::getTimeEnd()
{
    return this->timeEnd;
}

double NumericalSolveOfMotionEquationBase::getStepLength()
{
    return this->stepLength;
}

NumericalSolveOfMotionEquationBase::~NumericalSolveOfMotionEquationBase()
{
    this->position.clear();
    this->velocity.clear();
    this->momentum.clear();
    this->energy.clear();
    this->time.clear();
}