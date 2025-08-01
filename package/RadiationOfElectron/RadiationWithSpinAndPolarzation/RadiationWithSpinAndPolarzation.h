#include "../BasicRadiation/BasicRadiation.h"
#include "../../Dimension3Vector/Dimension3Vector.h"
#include <complex> 

class RadiationWithSpinAndPolarzation : public BasicRadiationOfElectronInCounterpropagatingLaser {
    private:
        double spinIncident;
        double spinEmission;
        double polarizationAlpha;
        double polarizationBeta;
        double emissionPolarAngleMin;
        double emissionPolarAngleMax;
        double axisOfIncidentAzimuthalAngleOfElectronSpin;
        double axisOfIncidentPolarAngleOfElectronSpin;
        double axisOfEmissionAzimuthalAngleOfElectronSpin;
        double axisOfEmissionPolarAngleOfElectronSpin;
        Dimension3Vector<double> incidentOrientationAxis;
        Dimension3Vector<double> emissionOrientationAxis;
        Dimension3Vector<double> combinedIncidentOrientationAxis;
        Dimension3Vector<double> combinedEmissionOrientationAxis;
        std::vector<double> stokesParameter;
    public:
        RadiationWithSpinAndPolarzation(double momentumZPrime,double momentumXPrime,double fieldParameter1,double fieldParameter2,double properTime,double photonEnergy,double emissionAzimuthalAngle,double spinIncident,double spinEmission,double polarizationAlpha,double polarzationBeta,double axisOfIncidentAzimuthalAngleOfElectronSpin,double axisOfIncidentPolarAngleOfElectronSpin,double axisOfEmissionAzimuthalAngleOfElectronSpin,double axisOfEmissionPolarAngleOfElectronSpin, double rotationDirection1, double rotationDirection2,double emissionPolarAngleMin,double emissionPolarAngleMax);
        //the value of the spin are 0.5 or -0.5
        void calculateDifferentialEmissionIntensity();
        void calculateDifferentialEmissionIntensityWithDoubledLabel();
        void calculateVortexDifferentialEmissionIntensity(double angularQuantumNumber, double polarizationParameter, size_t azimuthalAngleDivisions);
        void calculateStokesParameter();
        std::vector<double> calculateSixTermsOfDifferentialEmissionIntensity();
        double getSpinIncident();
        double getSpinEmission();
        double getAxisOfIncidentAzimuthalAngleOfElectronSpin();
        double getAxisOfIncidentPolarAngleOfElectronSpin();
        double getAxisOfEmissionAzimuthalAngleOfElectronSpin();
        double getAxisOfEmissionPolarAngleOfElectronSpin();
        Dimension3Vector<double> getIncidentOrientationAxis();
        Dimension3Vector<double> getEmissionOrientationAxis();
        std::vector<double> getStokesParameter();
        std::vector<double> getStokesParameterNormalized();
        ~RadiationWithSpinAndPolarzation();
};