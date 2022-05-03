/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_OCEAN_HPP
#define CPP_OCEAN_HPP

#include "Globals.hpp"
#include "Pirf.hpp"
#include "Prop.hpp"

#include <iostream>
#include <string>

using namespace std;

class Ocean {
public:
    Ocean();

    Ocean(Pirf pirf, double hmix, double cp, double dens, double dens_c, double aoc, double kg_aoc, double t_chem,
          string_view model_name, string doc);

    Ocean(Ocean &&ocean) noexcept;

    Ocean &operator=(Ocean &&ocean) noexcept;

    ~Ocean();

    void initCarbonStock(int ntime);

    void increaseTempKElem(int i, double value);

    void setDpCsElem(int i, double value);

    void setMsk(double *new_array);

    void increaseMskElem(int i, double value);

    void setMsElem(int i, double value);

    void setTempk(double *new_array);

    void resetTempk();

    void setTempElem(int i, double value);

    void co2Exceed();

    void tExceed();

    [[nodiscard]] double getAoc() const;

    [[nodiscard]] double getKgAoc() const;

    [[nodiscard]] double getOmt() const;

    [[nodiscard]] double getOmc() const;

    [[nodiscard]] string getDoc() const;

    [[nodiscard]] double getDpCsElem(int i) const;

    Prop &getPropagator();

    [[nodiscard]] int getPropNscale() const;

    [[nodiscard]] double getPropfElem(int i) const;

    double *getMsk();

    [[nodiscard]] double getMskElem(int i) const;

    [[nodiscard]] double getMsElem(int i) const;

    double *getTempK();

    [[nodiscard]] double getTempKElem(int i) const;

    [[nodiscard]] double getTempElem(int i) const;

    [[nodiscard]] bool isCo2Exceed() const;

    [[nodiscard]] bool isTExceed() const;

    void computePropagators(double dt);

    double computeDpCo2s(bool t_dep, double co2_atm0, int i, double t, bool deriv);

private:
    Pirf _pirf;
    double _hmix = 0e0;
    double _cp = 0e0;
    double _dens = 0e0;
    double _dens_c = 0e0;
    double _aoc = 0e0;
    double _kg_aoc = 0e0;
    double _t_chem = 0e0;
    double _om_t = 0e0;
    double _om_c = 0e0;
    string _model_name;
    string _doc;
    double const _buffer_t = 0.0423e0;
    double *_dpcs{};
    Prop _prop = Prop(0, 0e0);
    double *_msk{};
    double *_ms{};
    double *_temp{};
    double *_temp_k{};
    bool _co2_exceed = false;
    bool _t_exceed = false;
};


#endif //CPP_OCEAN_HPP
