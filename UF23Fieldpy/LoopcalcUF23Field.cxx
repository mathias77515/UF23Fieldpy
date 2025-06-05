#include "LoopcalcUF23Field.h"
#include "UF23Field.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void LoopUF23Field::printProgressBar(int progress, int total, bool display) {
    const int barWidth = 50;
    if (display) {
        float percentage = static_cast<float>(progress) / total * 100;

        std::cout << "[";
        int pos = barWidth * progress / total;
        for (int i = 0; i < barWidth; ++i)
            std::cout << (i < pos ? "#" : " ");
        std::cout << "] " << int(percentage) << " %\r";
        std::cout.flush();
    }
    
}

Vector3 LoopUF23Field::sphericalToCartesian(double r, double theta, double phi) {
    return Vector3(r * sin(theta) * cos(phi),
                   r * sin(theta) * sin(phi),
                   r * cos(theta));
}

Vector3 LoopUF23Field::cartesianFieldToSpherical(double Bx, double By, double Bz, double theta, double phi) {
    double sin_theta = sin(theta), cos_theta = cos(theta);
    double sin_phi = sin(phi), cos_phi = cos(phi);
    double Br     = sin_theta * cos_phi * Bx + sin_theta * sin_phi * By + cos_theta * Bz;
    double Btheta = cos_theta * cos_phi * Bx + cos_theta * sin_phi * By - sin_theta * Bz;
    double Bphi   = -sin_phi * Bx + cos_phi * By;
    return Vector3(Br, Btheta, Bphi);
}

LoopUF23Field::LoopUF23Field(int nside, int nlayers, const std::vector<double>& rValues, const std::vector<double>& thetaValues, const std::vector<double>& phiValues, UF23Field::ModelType model_type)
    : nside_(nside), nlayers_(nlayers), rValues_(rValues), thetaValues_(thetaValues), phiValues_(phiValues),
      uf23Field_(model_type) {
    npix_ = 12 * nside_ * nside_;
    computeDeltaR();
}

void LoopUF23Field::set_parameters(const std::vector<double>& params) {
    uf23Field_.SetParameters(params);
}
std::vector<double> LoopUF23Field::get_parameters() {
    return uf23Field_.GetParameters();

}
void LoopUF23Field::loadThetaPhiFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) throw std::runtime_error("Error opening file: " + filename);

    std::string line;
    while (std::getline(inFile, line)) {
        std::istringstream ss(line);
        double theta, phi;
        if (ss >> theta >> phi) {
            thetaValues_.push_back(theta);
            phiValues_.push_back(phi);
        }
    }
}

void LoopUF23Field::computeDeltaR() {
    delta_rValues_.resize(nlayers_ - 1);
    for (int i = 0; i < nlayers_ - 1; ++i) {
        delta_rValues_[i] = (rValues_[i+1] - rValues_[i]) / (nlayers_ - 1);
    }
    delta_rValues_.push_back(0.);  // add last layer step size as 0 (optional)
}

std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>
> LoopUF23Field::operator()(double offset_x=0, double offset_y=0, double offset_z=0, bool display=true) {

    std::vector<std::vector<double>> radial_distance_value(nlayers_, std::vector<double>(npix_, 0.0));
    std::vector<std::vector<double>> theta_value(nlayers_, std::vector<double>(npix_, 0.0));
    std::vector<std::vector<double>> phi_value(nlayers_, std::vector<double>(npix_, 0.0));

    std::vector<std::vector<double>> Br_value(nlayers_, std::vector<double>(npix_, 0.0));
    std::vector<std::vector<double>> Btheta_value(nlayers_, std::vector<double>(npix_, 0.0));
    std::vector<std::vector<double>> Bphi_value(nlayers_, std::vector<double>(npix_, 0.0));
    
    for (int pixel = 0; pixel < npix_; ++pixel) {
        for (int i = 0; i < nlayers_; ++i) {
            Vector3 coords_cart = sphericalToCartesian(rValues_[i], thetaValues_[pixel], phiValues_[pixel]);
            coords_cart.x += offset_x;
            coords_cart.y += offset_y;
            coords_cart.z += offset_z;

            const auto field = uf23Field_(coords_cart);
            const auto field_sph = cartesianFieldToSpherical(field.x, field.y, field.z, thetaValues_[pixel], phiValues_[pixel]);

            radial_distance_value[i][pixel] = rValues_[i];
            theta_value[i][pixel] = thetaValues_[pixel];
            phi_value[i][pixel] = phiValues_[pixel];

            Br_value[i][pixel] = field_sph.x;
            Btheta_value[i][pixel] = field_sph.y;
            Bphi_value[i][pixel] = field_sph.z;
            
        }
        
        // Ensures only one thread writes to stdout at a time (avoids jumbled output)
        printProgressBar(pixel, npix_, display);
    }

    return std::make_tuple(
        std::move(radial_distance_value),
        std::move(theta_value),
        std::move(phi_value),
        std::move(Br_value),
        std::move(Btheta_value),
        std::move(Bphi_value)
    );
}


// Bind Vector3 to Python
PYBIND11_MODULE(UF23Fieldpy, m) {
  py::class_<Vector3>(m, "Vector3")
    //.def(py::init<>())
    .def(py::init<double, double, double>())
    .def_readwrite("x", &Vector3::x)
    .def_readwrite("y", &Vector3::y)
    .def_readwrite("z", &Vector3::z);

  // Bind UF23Field::ModelType enum (assuming it exists as enum ModelType)
  py::enum_<UF23Field::ModelType>(m, "ModelType")
    .value("base", UF23Field::ModelType::base)
    .value("cre10", UF23Field::ModelType::cre10)
    .value("nebCor", UF23Field::ModelType::nebCor)
    .value("neCL", UF23Field::ModelType::neCL)
    .value("spur", UF23Field::ModelType::spur)
    .value("synCG", UF23Field::ModelType::synCG)
    .value("twistX", UF23Field::ModelType::twistX)
    .value("expX", UF23Field::ModelType::expX)
    .export_values();

  py::class_<UF23Field>(m, "UF23Field")
    .def(py::init<UF23Field::ModelType, double>(), py::arg("model_type"), py::arg("max_radius_kpc"))
    // Bind operator() as a method named "evaluate" (or __call__ for Python callable)
    .def("get_parameters", &UF23Field::GetParameters)
    .def("set_parameters", &UF23Field::SetParameters)
    .def("__call__", [](const UF23Field &self, const Vector3 &pos) {
      return self(pos);
    }, py::arg("position"));

  // Add this for LoopUF23Field:
  py::class_<LoopUF23Field>(m, "LoopUF23Field")
    .def(py::init<int, int, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, UF23Field::ModelType>(), 
         py::arg("nside"), py::arg("nlayers"), py::arg("rValues"), py::arg("thetaValues"), py::arg("phiValues"), py::arg("model_type") = UF23Field::ModelType::base)
    .def("set_parameters", &LoopUF23Field::set_parameters, py::arg("params"))
    .def("get_parameters", &LoopUF23Field::get_parameters)
    .def("__call__", &LoopUF23Field::operator());
}