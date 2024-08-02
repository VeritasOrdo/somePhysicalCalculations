#include<vector>
#include<cmath>
#include<complex>   

#pragma once

template <typename T>

class Dimension3Vector {
    private:
        std::vector<T> vector;
    public:
        Dimension3Vector();
        Dimension3Vector(T x, T y, T z);
        Dimension3Vector(const Dimension3Vector<T>& vector);
        ~Dimension3Vector();
        T getX() const;
        T getY() const;
        T getZ() const;
        void setX(T x);
        void setY(T y);
        void setZ(T z);
        //get the length of the vector(complex number)
        long double getLength() const;
        Dimension3Vector<T> operator+(const Dimension3Vector<T>& vector) const;
        Dimension3Vector<T> operator-(const Dimension3Vector<T>& vector) const;
        Dimension3Vector<T> operator*(T scalar) const;
        Dimension3Vector<T> operator/(T scalar) const;
        T operator*(const Dimension3Vector<T>& vector) const;
        //std::complex<double> operator*(const Dimension3Vector<double>& vector) const;
        Dimension3Vector<T> operator^(const Dimension3Vector<T>& vector) const;
        Dimension3Vector<T> operator-() const;
        bool operator==(const Dimension3Vector<T>& vector) const;
        bool operator!=(const Dimension3Vector<T>& vector) const;
        Dimension3Vector<T>& operator=(const Dimension3Vector<T>& vector);
        Dimension3Vector<T>& operator+=(const Dimension3Vector<T>& vector);
        Dimension3Vector<T>& operator-=(const Dimension3Vector<T>& vector);
        Dimension3Vector<T>& operator*=(T scalar);
        Dimension3Vector<T>& operator/=(T scalar);
        T& operator[](int index);
        const T& operator[](int index) const;
        void normalize();
        Dimension3Vector<T> normalized() const;
        Dimension3Vector<T> conjugate() const;
        operator Dimension3Vector<std::complex<double>>() const;
        friend std::ostream& operator<<(std::ostream& os, const Dimension3Vector<T>& vector) {
            os << "(" << vector.vector[0] << "," << vector.vector[1] << "," << vector.vector[2] << ")";
            return os;
        }
};