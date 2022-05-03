/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_ATMOSPHERE_HPP
#define CPP_ATMOSPHERE_HPP

#include "Globals.hpp"
#include "UtilsMath.hpp"

class Atmosphere {
public:
    Atmosphere();

    Atmosphere(Atmosphere &&atmosphere) noexcept;

    Atmosphere &operator=(Atmosphere &&atmosphere) noexcept;

    ~Atmosphere();

    void initialize(int ntime);

    void initializeCo2Atm0();

    void setMaElem(int i, double value);

    void setECo2Elem(int i, double value);

    double getMaElem(int i) const;

    double getECo2Elem(int i) const;

    double getCo2Atm0() const;

    double rfCo2(int i);

private:
    double *_ma{};
    double *_e_co2{};
    double _co2_atm0 = 0e0;
};


#endif //CPP_ATMOSPHERE_HPP
