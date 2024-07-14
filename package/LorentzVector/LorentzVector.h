#include<vector>
#include<string>
#include<complex>

class LorentzVector {
    private:
        double components[4];
    public:
        LorentzVector(double t, double x, double y, double z);
        double get(int i);
        double &operator[](int i);
        double dot(LorentzVector other);
        double operator*(LorentzVector other);
        LorentzVector add(LorentzVector other);
        LorentzVector operator+(LorentzVector other);
        LorentzVector multiply(double scalar);
        std::complex<double> magnitude();
        std::string toString();
};

