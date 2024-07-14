#include "LorentzVector.h"

LorentzVector::LorentzVector(double t, double x, double y, double z) {
    components[0] = t;
    components[1] = x;
    components[2] = y;
    components[3] = z;
}

double LorentzVector::get(int i) {
    return components[i];
}

double &LorentzVector::operator[](int i) {
    return components[i];
}

double LorentzVector::dot(LorentzVector other) {
    return components[0]*other[0] - components[1]*other[1] - components[2]*other[2] - components[3]*other[3];
}

double LorentzVector::operator*(LorentzVector other) {
    return dot(other);
}

LorentzVector LorentzVector::add(LorentzVector other) {
    return LorentzVector(components[0] + other[0], components[1] + other[1], components[2] + other[2], components[3] + other[3]);
}

LorentzVector LorentzVector::operator+(LorentzVector other) {
    return add(other);
}

LorentzVector LorentzVector::multiply(double scalar) {
    return LorentzVector(components[0]*scalar, components[1]*scalar, components[2]*scalar, components[3]*scalar);
}

std::complex<double> LorentzVector::magnitude() {
    return std::sqrt(components[0]*components[0] - components[1]*components[1] - components[2]*components[2] - components[3]*components[3]);
}

std::string LorentzVector::toString() {
    return "(" + std::to_string(components[0]) + ", " + std::to_string(components[1]) + ", " + std::to_string(components[2]) + ", " + std::to_string(components[3]) + ")";
}