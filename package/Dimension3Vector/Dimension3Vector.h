#include<vector>
#include<cmath>

class Dimension3Vector{
    private:
        std::vector<double> vector;
    public:
        Dimension3Vector();
        Dimension3Vector(double x,double y,double z);
        void setVector(double x,double y,double z);
        void setVector(std::vector<double> vector);
        std::vector<double> getVector();
        double getNorm();
        double getNormSquare();
        double getDotProduct(Dimension3Vector vector);
        Dimension3Vector getCrossProduct(Dimension3Vector vector);
        Dimension3Vector getUnitVector();
        Dimension3Vector operator+(Dimension3Vector vector);
        Dimension3Vector operator-(Dimension3Vector vector);
        Dimension3Vector operator*(double scalar);
        Dimension3Vector operator/(double scalar);
        Dimension3Vector operator+=(Dimension3Vector vector);
        Dimension3Vector operator-=(Dimension3Vector vector);
        Dimension3Vector operator*=(double scalar);
        Dimension3Vector operator/=(double scalar);
        bool operator==(Dimension3Vector vector);
        bool operator!=(Dimension3Vector vector);
        ~Dimension3Vector();
};