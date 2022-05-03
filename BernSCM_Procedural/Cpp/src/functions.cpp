//
// Created by Brice Petit on 25-02-22.
//

#include "functions.hpp"

/**
 * Do the interpolation. It is the process of estimating unknown values that
 * fall between known values.
 *
 * @param x     X-value for interpolation.
 * @param x0    X-value for interpolation.
 * @param x1    X-value for interpolation.
 * @param y0    Y(x) interpol value.
 * @param y1    Y(x) interpol value.
 * @return      Return the interpol value.
 */
double interpol(double &x, double &x0, double &x1, double &y0, double &y1) {
    double y;
    if (y0 == NA) {
        if (abs(x - x1) < 1e-9) {
            y = y1;
        } else {
            y = NA;
        }
    } else if (y1 == NA) {
        if (abs(x - x0) < 1e-9) {
            y = y0;
        } else {
            y = NA;
        }
    } else {
        y = ((x - x0) * y1 + (x1 - x) * y0) / (x1 - x0);
    }
    return y;
}

/**
 * Do the computation of the npp according to the model.
 *
 * @param _ma   Mass of C in atmosphere.
 * @param t     Global ΔSAT (℃)
 * @param deriv Derivative dNPP/dm.
 * @return      Return the npp.
 */
double npp(double _ma, double t, bool deriv) {
    double res = 0e0;
    if (LAND_MODEL == "HRBM") {
        res = computeHrbmNpp(_ma, t, deriv);
    } else if (LAND_MODEL == "4box") {
        res = compute4boxNpp(_ma, deriv);
    }
    return res;
}
/**
 * Compute the ne
 *
 * t primary production for the model HRBM.
 *
 * @param _ma       Mass of C in atmosphere.
 * @param t         Global ΔSAT (℃)
 * @param deriv     Derivative dNPP/dm.
 * @return          Return the npp.
 */
double computeHrbmNpp(double _ma, double t, bool deriv) {
    double npp;
    if (_ma / PPMTOGT > 1120 && !co2_exceed) {
        cerr << "warning: CO2 parametrization range for NPP in HRBM exceeded (1120ppm)" << endl;
        co2_exceed = true;
    }
    if (t > 5 && !t_exceed) {
        cerr << "warning: temperature parametrization range for NPP in HRBM exceeded (5K)" << endl;
        t_exceed = true;
    }
    if (!deriv) {
        if (co2_dep) {
            npp = -exp(3.672801e0) + exp(-0.430818e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 1) -
                  exp(-6.145559e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 2) +
                  exp(-12.353878e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 3) -
                  exp(-19.010800e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 4) +
                  exp(-26.183752e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 5) -
                  exp(-34.317488e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 6) -
                  exp(-41.553715e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 7) +
                  exp(-48.265138e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 8) -
                  exp(-56.056095e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 9) +
                  exp(-64.818185e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 10);
        } else {
            npp = NPP0;
        }
        if (t_dep) {
            npp *= (1 + 0.11780208e+0 * tanh(t / 0.509312421e+02) + 0.24305130e-02 * tanh(t / 0.885326739e+01));
        }
    } else {
        double npp_dev;
        if (co2_dep) {
            npp_dev = (+exp(-0.430818e0) - 2 * exp(-6.145559e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 1) +
                       3 * exp(-12.353878e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 2) -
                       4 * exp(-19.010800e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 3) +
                       5 * exp(-26.183752e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 4) -
                       6 * exp(-34.317488e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 5) -
                       7 * exp(-41.553715e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 6) +
                       8 * exp(-48.265138e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 7) -
                       9 * exp(-56.056095e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 8) +
                       10 * exp(-64.818185e0 + myLog(min(_ma / PPMTOGT, 1274e0)) * 9)) / PPMTOGT;
        } else {
            npp_dev = 0e0;
        }
        if (t_dep) {
            npp_dev *= (1e0 + 0.11780208e+0 * tanh(t / 0.509312421e+02) + 0.24305130e-02 * tanh(t / 0.885326739e+01));
        }
        npp = npp_dev;
    }
    return npp;
}

/**
 * Compute the net primary production for the 4box model.
 *
 * @param _ma       Mass of C in atmosphere.
 * @param deriv     Derivative dNPP/dm.
 * @return          Return the npp.
 */
double compute4boxNpp(double _ma, bool deriv) {
    double npp, m0, d_npp;
    if (!deriv) {
        if (co2_dep) {
            m0 = co2_atm0 * PPMTOGT;
            d_npp = NPP0 * FERT * myLog(_ma / m0);
            npp = d_npp + NPP0;
        } else {
            npp = NPP0;
        }
    } else {
        double npp_dev;
        if (co2_dep) {
            npp_dev = NPP0 * FERT * (1 / _ma);
        } else {
            npp_dev = 0e0;
        }
        npp = npp_dev;
    }
    return npp;
}

/**
 * We create our specific log because the owner use a certain flag in order to do the following.
 *
 * @param value The value to use in the logarithm.
 * @return      Return the computed value.
 */
double myLog(const double &value) {
    double res;
    if (value > 0e0) {
        res = log(value);
    } else {
        res = value;
    }
    return res;
}

/**
 *
 * @param ca        Atmospheric CO2 concentration (Gt).
 * @param dp_co2s   Ocean saturation CO2 pressure deviation from preindustrial equilibrium (ppm).
 * @return          Return the atmosphere-ocean CO2 flux (Gt/yr)
 */
double fasC(double &ca, double &dp_co2s) {
    return kg_aoc * ((ca - co2_atm0 * PPMTOGT) - dp_co2s * PPMTOGT);
}

/**
 * Compute the air-sea heat flux.
 *
 * @param rf_tot    Radiative forcing (Wm⁻²).
 * @param t         Global near-surface atmospheric temperature deviation (℃).
 * @return          Return the air-sea heat flux (PW).
 */
double fasT(double &rf_tot, double &t) {
    if (t2x > 0e0) {
        return (aoc / OFRAC / PETA) * (rf_tot - (t / t2x) * RF2X);
    } else {
        // T is always in equilibrium
        return 0e0;
    }
}

/**
 * Compute the RF of the atmospheric CO2.
 *
 * @param _ma   Atmospheric CO₂ (Gt).
 * @return      Return RF of atmospheric CO2 (Wm⁻²).
 */
double rfCo2(double &_ma) {
    return RECO2 * myLog((_ma / PPMTOGT) / CO2PREIND);
}

/**
 * Calculate equivalent atmospheric CO2 (in GtC) from RF.
 *
 * @param rf_co2    Co2 RF (Wm⁻²).
 * @return          Return the atmospheric CO₂ (Gt).
 */
double rfEqC02Ma(double rf_co2) {
    return exp(rf_co2 / RECO2) * CO2PREIND * PPMTOGT;
}

/**
 * Analytical representation of the zeta-factor for a temperature
 * range from 17.7 to 18.3 degrees celsius (Tchem). Rogers chemistry model
 * was used to calculate the zeta-factor.
 * @param d_dic Change in ocean surface DIC (μmol:kg).
 * @param t     Global SAT change from preindustrial (℃).
 * @param deriv Derivative dpCO2s/ddDIC.
 * @return      Ocean saturation CO2 pressure deviation from preindustrial
 *              equilibrium (ppm), or derivative (d dpCs/d dDIC).
 */
double computeDpCo2s(double d_dic, double &t, bool deriv) {
    double dp_co2s;
    if (!deriv) {
        // Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
        dp_co2s = (1.5568e0 - 1.3993e-2 * t_chem) * d_dic +
                  (7.4706e0 - 0.20207e0 * t_chem) * 1e-3 * pow(d_dic, 2) -
                  (1.2748e0 - 0.12015e0 * t_chem) * 1e-5 * pow(d_dic, 3) +
                  (2.4491e0 - 0.12639e0 * t_chem) * 1e-7 * pow(d_dic, 4) -
                  (1.5468e0 - 0.15326e0 * t_chem) * 1e-10 * pow(d_dic, 5);
        if (dp_co2s > 1320 && !co2_exceed) {
            cerr << "warning: CO2 parametrization range for dpCO2s exceeded (1320ppm)" << endl;
            co2_exceed = true;
        }
        if (t + t_chem > 25 && !t_exceed) {
            cerr << "warning: temperature parametrization range for dpCO2s exceeded (25℃)" << endl;
            t_exceed = true;
        }
        if (t_dep) {
            dp_co2s = (dp_co2s + co2_atm0) * exp(buffer_t * t) - co2_atm0;
        }
    } else {
        double dp_co2s_dev;
        dp_co2s_dev =
                (1.5568e0 - 1.3993e-2 * t_chem) + 2e0 * (7.4706 - 0.20207e0 * t_chem) * 1.0e-3 * d_dic -
                3e0 * (1.2748 - 0.12015e0 * t_chem) * 1.0e-5 * pow(d_dic, 2) +
                4e0 * (2.4491 - 0.12639e0 * t_chem) * 1.0e-7 * pow(d_dic, 3) -
                5e0 * (1.5468 - 0.15326e0 * t_chem) * 1.0e-10 * pow(d_dic, 4);
        if (t_dep) {
            dp_co2s_dev = dp_co2s_dev * exp(buffer_t * t);
        }
        dp_co2s = dp_co2s_dev;
    }
    return dp_co2s;
}
