/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_SIRF_HPP
#define CPP_SIRF_HPP

#include "Globals.hpp"

#include <cmath>
#include <cstring>

using namespace std;

class Sirf {
public:
    Sirf();

    Sirf(double t_max, double weight[], double t_scale[]);

    Sirf(const Sirf &sirf);

    Sirf(Sirf &&sirf) noexcept;

    Sirf &operator=(Sirf &&sirf) noexcept;

    ~Sirf();

    [[nodiscard]] double getTMax() const;

    [[nodiscard]] double getWeightElem(int i) const;

    [[nodiscard]] double getTscaleElem(int i) const;

private:
    double _t_max = 0e0;
    double *_weight{};
    double *_t_scale{};
};


#endif //CPP_SIRF_HPP
