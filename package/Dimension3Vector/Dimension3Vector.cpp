#include "Dimension3Vector.h"

template <>
Dimension3Vector<double>::operator Dimension3Vector<std::complex<double>>() const {
    return Dimension3Vector<std::complex<double>>(this->getX(),this->getY(),this->getZ());
}

template <>
Dimension3Vector<int> Dimension3Vector<int>::conjugate() const {
    return *this;
}

template<>
Dimension3Vector<float> Dimension3Vector<float>::conjugate() const {
    return *this;
}

template<>
Dimension3Vector<double> Dimension3Vector<double>::conjugate() const {
    return *this;
}

template<>
Dimension3Vector<long double> Dimension3Vector<long double>::conjugate() const {
    return *this;
}

template<>
Dimension3Vector<std::complex<int>> Dimension3Vector<std::complex<int>>::conjugate() const {
    return Dimension3Vector<std::complex<int>>(std::conj(this->getX()),std::conj(this->getY()),std::conj(this->getZ()));
}

template<>
Dimension3Vector<std::complex<float>> Dimension3Vector<std::complex<float>>::conjugate() const {
    return Dimension3Vector<std::complex<float>>(std::conj(this->getX()),std::conj(this->getY()),std::conj(this->getZ()));
}

template<>
Dimension3Vector<std::complex<double>> Dimension3Vector<std::complex<double>>::conjugate() const {
    return Dimension3Vector<std::complex<double>>(std::conj(this->getX()),std::conj(this->getY()),std::conj(this->getZ()));
}

template<>
Dimension3Vector<std::complex<long double>> Dimension3Vector<std::complex<long double>>::conjugate() const {
    return Dimension3Vector<std::complex<long double>>(std::conj(this->getX()),std::conj(this->getY()),std::conj(this->getZ()));
}

template <typename T>
Dimension3Vector<T>::Dimension3Vector() {
    vector.push_back(0);
    vector.push_back(0);
    vector.push_back(0);
}

template <typename T>
Dimension3Vector<T>::Dimension3Vector(T x, T y, T z) {
    vector.push_back(x);
    vector.push_back(y);
    vector.push_back(z);
}

template <typename T>
Dimension3Vector<T>::Dimension3Vector(const Dimension3Vector<T>& vector) {
    this->vector = vector.vector;
}

template <typename T>
Dimension3Vector<T>::~Dimension3Vector() {
    vector.clear();
}

template <typename T>
T Dimension3Vector<T>::getX() const {
    return vector[0];
}

template <typename T>
T Dimension3Vector<T>::getY() const {
    return vector[1];
}

template <typename T>
T Dimension3Vector<T>::getZ() const {
    return vector[2];
}

template <typename T>
void Dimension3Vector<T>::setX(T x) {
    vector[0] = x;
}

template <typename T>
void Dimension3Vector<T>::setY(T y) {
    vector[1] = y;
}

template <typename T>
void Dimension3Vector<T>::setZ(T z) {
    vector[2] = z;
}

template <typename T>
long double Dimension3Vector<T>::getLength() const {
    return sqrt(std::abs(vector[0]) * std::abs(vector[0]) + std::abs(vector[1]) * std::abs(vector[1]) + std::abs(vector[2]) * std::abs(vector[2]));
}

template <typename T>
Dimension3Vector<T> Dimension3Vector<T>::operator+(const Dimension3Vector<T>& vector) const {
    return Dimension3Vector<T>(this->vector[0] + vector.vector[0], this->vector[1] + vector.vector[1], this->vector[2] + vector.vector[2]);
}

template <typename T>
Dimension3Vector<T> Dimension3Vector<T>::operator-(const Dimension3Vector<T>& vector) const {
    return Dimension3Vector<T>(this->vector[0] - vector.vector[0], this->vector[1] - vector.vector[1], this->vector[2] - vector.vector[2]);
}

template<typename T>
Dimension3Vector<T> Dimension3Vector<T>::operator*(T scalar) const {
    return Dimension3Vector<T>(this->vector[0] * scalar, this->vector[1] * scalar, this->vector[2] * scalar);
}

template <typename T>
Dimension3Vector<T> Dimension3Vector<T>::operator/(T scalar) const {
    return Dimension3Vector<T>(this->vector[0] / scalar, this->vector[1] / scalar, this->vector[2] / scalar);
}

template <typename T>
T Dimension3Vector<T>::operator*(const Dimension3Vector<T>& vector) const {
    return this->vector[0] * vector.vector[0] + this->vector[1] * vector.vector[1] + this->vector[2] * vector.vector[2];
}

template <typename T>
Dimension3Vector<T> Dimension3Vector<T>::operator^(const Dimension3Vector<T>& vector) const {
    return Dimension3Vector<T>(this->vector[1] * vector.vector[2] - this->vector[2] * vector.vector[1], this->vector[2] * vector.vector[0] - this->vector[0] * vector.vector[2], this->vector[0] * vector.vector[1] - this->vector[1] * vector.vector[0]);
}

template <typename T>
Dimension3Vector<T> Dimension3Vector<T>::operator-() const {
    return Dimension3Vector<T>(-this->vector[0], -this->vector[1], -this->vector[2]);
}

template <typename T>
bool Dimension3Vector<T>::operator==(const Dimension3Vector<T>& vector) const {
    return this->vector[0] == vector.vector[0] && this->vector[1] == vector.vector[1] && this->vector[2] == vector.vector[2];
}

template <typename T>
bool Dimension3Vector<T>::operator!=(const Dimension3Vector<T>& vector) const {
    return this->vector[0] != vector.vector[0] || this->vector[1] != vector.vector[1] || this->vector[2] != vector.vector[2];
}

template <typename T>
Dimension3Vector<T>& Dimension3Vector<T>::operator=(const Dimension3Vector<T>& vector) {
    this->vector = vector.vector;
    return *this;
}

template <typename T>
Dimension3Vector<T>& Dimension3Vector<T>::operator+=(const Dimension3Vector<T>& vector) {
    this->vector[0] += vector.vector[0];
    this->vector[1] += vector.vector[1];
    this->vector[2] += vector.vector[2];
    return *this;
}

template <typename T>
Dimension3Vector<T>& Dimension3Vector<T>::operator-=(const Dimension3Vector<T>& vector) {
    this->vector[0] -= vector.vector[0];
    this->vector[1] -= vector.vector[1];
    this->vector[2] -= vector.vector[2];
    return *this;
}

template <typename T>
Dimension3Vector<T>& Dimension3Vector<T>::operator*=(T scalar) {
    this->vector[0] *= scalar;
    this->vector[1] *= scalar;
    this->vector[2] *= scalar;
    return *this;
}

template <typename T>
Dimension3Vector<T>& Dimension3Vector<T>::operator/=(T scalar) {
    this->vector[0] /= scalar;
    this->vector[1] /= scalar;
    this->vector[2] /= scalar;
    return *this;
}

template <typename T>
T& Dimension3Vector<T>::operator[](int index) {
    return vector[index];
}

template <typename T>
const T& Dimension3Vector<T>::operator[](int index) const {
    return vector[index];
}

template <typename T>
void Dimension3Vector<T>::normalize() {
    T length = getLength();
    vector[0] /= length;
    vector[1] /= length;
    vector[2] /= length;
}

template <typename T>
Dimension3Vector<T> Dimension3Vector<T>::normalized() const {
    T length = getLength();
    return Dimension3Vector<T>(vector[0] / length, vector[1] / length, vector[2] / length);
}

template class Dimension3Vector<int>;
template class Dimension3Vector<float>;
template class Dimension3Vector<double>;
template class Dimension3Vector<long double>;
template class Dimension3Vector<std::complex<int>>;
template class Dimension3Vector<std::complex<float>>;
template class Dimension3Vector<std::complex<double>>;
template class Dimension3Vector<std::complex<long double>>;
