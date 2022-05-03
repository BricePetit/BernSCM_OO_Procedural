/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_EARTH_HPP
#define CPP_EARTH_HPP

#include "Atmosphere.hpp"
#include "Globals.hpp"
#include "Land.hpp"
#include "Ocean.hpp"
#include "Prop.hpp"

#include <string>

using namespace std;

class Earth {
public:
    Earth();

    Earth(Ocean &ocean, Land &land, Atmosphere &atmosphere, double **forcing, double t2x, bool t_dep, bool co2_dep,
          string scenario, double dt);

    Earth &operator=(Earth &&earth) noexcept;

    ~Earth();

    void initialize();

    void initAttributesToZero();

    void setForcing(int n);

    double interpol(int n, int index_forcing);

    void setRfbElem(int i, double value);

    void setFbElem(int i, double value);

    [[nodiscard]] string getOceanDoc() const;

    [[nodiscard]] double getTempElem(int i) const;

    [[nodiscard]] double getMsElem(int i) const;

    [[nodiscard]] string getLandDoc() const;

    [[nodiscard]] double getNpp0() const;

    [[nodiscard]] double getMlElem(int i) const;

    [[nodiscard]] double getMaElem(int i) const;

    [[nodiscard]] double getCo2Atm0() const;

    [[nodiscard]] double getT2x() const;

    [[nodiscard]] bool isTDep() const;

    [[nodiscard]] bool isCo2Dep() const;

    [[nodiscard]] string getScenario() const;

    [[nodiscard]] double getDt() const;

    [[nodiscard]] int getNtime() const;

    [[nodiscard]] double getTimeElem(int i) const;

    [[nodiscard]] double getRfElem(int i) const;

    [[nodiscard]] double getRfncElem(int i) const;

    [[nodiscard]] double getRfcElem(int i) const;

    [[nodiscard]] double getRfbElem(int i) const;

    [[nodiscard]] double getFhElem(int i) const;

    [[nodiscard]] double getFbElem(int i) const;

    [[nodiscard]] double getDpcsElem(int i) const;

    [[nodiscard]] double getFoElem(int i) const;

    [[nodiscard]] double *getFnpp() const;

    [[nodiscard]] double getFnppElem(int i) const;

    [[nodiscard]] double getECo2Elem(int i) const;

    void carbonCycleClimateSimulation();

    void timeStep(int n);

    double stepPulse(int n, const double *f, double *mk, Prop &q, double x);

    double fasC(int i);

    double fasT(int i);

    double npp(double ma, double t, bool deriv);

    double computeDpCo2s(int i, double t, bool deriv);

    template<typename S, typename T>
    void updateExceedValues(S &model, T &model_to_update);

private:
    Ocean _ocean;
    Land _land;
    Atmosphere _atmosphere;
    double **_forcing{};
    double _t2x = 0e0;
    bool _t_dep = false;
    bool _co2_dep = false;
    string _scenario;
    double _dt = 0e0;
    // Time dimension
    int _ntime = 0;
    // Time model
    double *_time{};
    // Forcing time series
    long int *_itime{};
    // Total radiative forcing (Wm⁻²)
    double *_rf{};
    // Non-CO₂ radiative forcing (Wm⁻²)
    double *_rfnc{};
    // CO₂ radiative forcing (Wm⁻²)
    double *_rfc{};
    // Budget radiative forcing (Wm⁻²)
    double *_rfb{};
    // Air-sea heat flux (PW)
    double *_fh{};
    // Budget C uptake (GtC/yr)
    double *_fb{};
    // Air-sea C flux (GtC/yr)
    double *_fo{};
    // NPP (GtC/yr)
    double *_fnpp{};
    // Budget C uptake boolean
    bool _f_budget = false;
    // Budget radiative forcing boolean
    bool _rf_budget = false;
    // CO2 budget boolean
    bool _co2_budget = false;
    bool _exceed = false;
    bool _t_exceed = false;
    bool _co2_exceed = false;
};


#endif //CPP_EARTH_HPP
