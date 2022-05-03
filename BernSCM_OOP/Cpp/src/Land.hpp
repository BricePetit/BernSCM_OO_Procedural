/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_LAND_HPP
#define CPP_LAND_HPP

#include "Globals.hpp"
#include "Pirf.hpp"
#include "Prop.hpp"
#include "Sirf.hpp"
#include "UtilsMath.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace std;

class Land {
public:
    Land();

    Land(Pirf pirf, Sirf sirf, double npp0, double fert, string_view model_name, string doc);

    Land(Pirf pirf, double npp0, double fert, string_view model_name, string doc);

    Land(Land &&Land) noexcept;

    Land &operator=(Land &&land) noexcept;

    ~Land();

    void initCarbonStock(int ntime);

    void increaseMlkElem(int index, double value);

    void setMlElem(int i, double value);

    void co2Exceed();

    void tExceed();

    [[nodiscard]] double getNpp0() const;

    [[nodiscard]] string getDoc() const;

    Prop &getPropagator();

    [[nodiscard]] int getPropNscale() const;

    [[nodiscard]] double getPropfElem(int i) const;

    double *getMlk();

    [[nodiscard]] double getMlElem(int i) const;

    [[nodiscard]] bool isExceed() const;

    [[nodiscard]] bool isCo2Exceed() const;

    [[nodiscard]] bool isTExceed() const;

    void computePropagators(double dt, Pirf *pirf = nullptr);

    void computePropTempDependent(double t, double dt);

    double computeNpp(double &ma, double &t, double const &co2_atm0, bool &co2_dep, bool &t_dep, bool &deriv);

    double computeHrbmNpp(double &ma, double &t, bool &co2_dep, bool &t_dep, bool &deriv);

    double compute4boxNpp(double &ma, double const &co2_atm0, bool &co2_dep, bool &deriv) const;

private:
    Pirf _pirf;
    Sirf _sirf;
    double _npp0 = 0e0;
    double _fert = 0e0;
    string _model_name;
    string _doc;
    Prop _prop = Prop(0, 0e0);
    double *_mlk{};
    double *_ml{};
    bool _exceed = false;
    bool _co2_exceed = false;
    bool _t_exceed = false;
};


#endif //CPP_LAND_HPP
