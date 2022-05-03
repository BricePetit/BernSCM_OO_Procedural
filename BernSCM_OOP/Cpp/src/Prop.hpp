/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_PROP_HPP
#define CPP_PROP_HPP

#include "Globals.hpp"
#include "Pirf.hpp"

#include <cmath>

class Prop {
public:
    Prop(int nscale, double x);

    Prop(Prop &&prop) noexcept;

    Prop &operator=(Prop &&prop) noexcept;

    ~Prop();

    [[nodiscard]] int getNscale() const;

    [[nodiscard]] double getPropmElem(int i) const;

    [[nodiscard]] double getPropfElem(int i) const;

    [[nodiscard]] double getPropfoElem(int i) const;

    void computePropagators(Pirf &pirf, double dt);

private:
    int _nscale;
    double _x;
    double _propm[NSCALEMAX + 1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    double _propf[NSCALEMAX + 1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    double _propfo[NSCALEMAX + 1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
};

#endif //CPP_PROP_HPP
