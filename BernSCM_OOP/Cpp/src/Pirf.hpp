/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_PIRF_HPP
#define CPP_PIRF_HPP

#include "Globals.hpp"

#include <cstring>
#include <string>

using namespace std;

class Pirf {
public:
    Pirf();

    Pirf(string_view name, int nscale, double weight[], double t_scale[]);

    Pirf(Pirf &&pirf) noexcept;

    Pirf(const Pirf &p) noexcept;

    Pirf &operator=(Pirf &&pirf) noexcept;

    ~Pirf();

    [[nodiscard]] string getName() const;

    [[nodiscard]] int getNscale() const;

    [[nodiscard]] double getWeightElem(int i) const;

    [[nodiscard]] double getTScaleElem(int i) const;

    void multiplyWeights(int bound, const double *value);

    void divideWeights(int bound, double value);

    void multiplyTscales(int bound, const double *value);

private:
    string _name;
    int _nscale = 0;
    double *_weight{};
    double *_t_scale{};
};


#endif //CPP_PIRF_HPP
