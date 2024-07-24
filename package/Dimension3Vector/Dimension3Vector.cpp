#include "Dimension3Vector.h"

Dimension3Vector::Dimension3Vector(){
    vector = {0,0,0};
}

Dimension3Vector::Dimension3Vector(double x,double y,double z){
    vector = {x,y,z};
}

void Dimension3Vector::setVector(double x,double y,double z){
    vector = {x,y,z};
}

void Dimension3Vector::setVector(std::vector<double> vector){
    this->vector = vector;
}

std::vector<double> Dimension3Vector::getVector(){
    return vector;
}

double Dimension3Vector::getNorm(){
    return std::sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
}

double Dimension3Vector::getNormSquare(){
    return vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2];
}

double Dimension3Vector::getDotProduct(Dimension3Vector vector){
    return this->vector[0]*vector.getVector()[0]+this->vector[1]*vector.getVector()[1]+this->vector[2]*vector.getVector()[2];
}

Dimension3Vector Dimension3Vector::getCrossProduct(Dimension3Vector vector){
    return Dimension3Vector(vector.getVector()[2]*this->vector[1]-vector.getVector()[1]*this->vector[2],vector.getVector()[0]*this->vector[2]-vector.getVector()[2]*this->vector[0],vector.getVector()[1]*this->vector[0]-vector.getVector()[0]*this->vector[1]);
}

Dimension3Vector Dimension3Vector::getUnitVector(){
    return *this/std::sqrt(this->getNormSquare());
}

Dimension3Vector Dimension3Vector::operator+(Dimension3Vector vector){
    return Dimension3Vector(this->vector[0]+vector.getVector()[0],this->vector[1]+vector.getVector()[1],this->vector[2]+vector.getVector()[2]);
}

Dimension3Vector Dimension3Vector::operator-(Dimension3Vector vector){
    return Dimension3Vector(this->vector[0]-vector.getVector()[0],this->vector[1]-vector.getVector()[1],this->vector[2]-vector.getVector()[2]);
}

Dimension3Vector Dimension3Vector::operator*(double scalar){
    return Dimension3Vector(this->vector[0]*scalar,this->vector[1]*scalar,this->vector[2]*scalar);
}

Dimension3Vector Dimension3Vector::operator/(double scalar){
    return Dimension3Vector(this->vector[0]/scalar,this->vector[1]/scalar,this->vector[2]/scalar);
}

Dimension3Vector Dimension3Vector::operator+=(Dimension3Vector vector){
    this->vector[0] += vector.getVector()[0];
    this->vector[1] += vector.getVector()[1];
    this->vector[2] += vector.getVector()[2];
    return *this;
}

Dimension3Vector Dimension3Vector::operator-=(Dimension3Vector vector){
    this->vector[0] -= vector.getVector()[0];
    this->vector[1] -= vector.getVector()[1];
    this->vector[2] -= vector.getVector()[2];
    return *this;
}

Dimension3Vector Dimension3Vector::operator*=(double scalar){
    this->vector[0] *= scalar;
    this->vector[1] *= scalar;
    this->vector[2] *= scalar;
    return *this;
}

Dimension3Vector Dimension3Vector::operator/=(double scalar){
    this->vector[0] /= scalar;
    this->vector[1] /= scalar;
    this->vector[2] /= scalar;
    return *this;
}

bool Dimension3Vector::operator==(Dimension3Vector vector){
    return this->vector[0] == vector.getVector()[0] && this->vector[1] == vector.getVector()[1] && this->vector[2] == vector.getVector()[2];
}

bool Dimension3Vector::operator!=(Dimension3Vector vector){
    return this->vector[0] != vector.getVector()[0] || this->vector[1] != vector.getVector()[1] || this->vector[2] != vector.getVector()[2];
}

Dimension3Vector::~Dimension3Vector(){
}

