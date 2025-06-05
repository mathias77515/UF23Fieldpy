#ifndef LOOPCALCUF23FIELD_H
#define LOOPCALCUF23FIELD_H

#include <vector>
#include <tuple>
#include "UF23Field.h"  // Make sure this includes Vector3 and UF23Field declarations

class LoopUF23Field {
public:
    LoopUF23Field(int nside, int nlayers, const std::vector<double>& rValues, const std::vector<double>& thetaValues, const std::vector<double>& phiValues, UF23Field::ModelType model_type = UF23Field::ModelType::base);

    // Change the operator() return type from void to a tuple of vectors
    std::tuple<
        std::vector<std::vector<double>>,
        std::vector<std::vector<double>>,
        std::vector<std::vector<double>>,
        std::vector<std::vector<double>>,
        std::vector<std::vector<double>>,
        std::vector<std::vector<double>>
    > operator()(double offset_x, double offset_y, double offset_z, bool display);
    void set_parameters(const std::vector<double>& params);
    std::vector<double> get_parameters();

private:
    int nside_;
    int nlayers_;
    int npix_;
    std::vector<double> rValues_;
    std::vector<double> thetaValues_;
    std::vector<double> phiValues_;
    std::vector<double> delta_rValues_;

    UF23Field uf23Field_;

    void loadThetaPhiFromFile(const std::string& filename);
    void computeDeltaR();

    // Helper functions for coordinate conversions
    Vector3 sphericalToCartesian(double r, double theta, double phi);
    Vector3 cartesianFieldToSpherical(double Bx, double By, double Bz, double theta, double phi);

    // Progress bar helper (optional)
    void printProgressBar(int progress, int total, bool display);
};

#endif // LOOPCALCUF23FIELD_H
