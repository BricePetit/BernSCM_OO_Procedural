//
// Created by Brice Petit on 25-02-22.
//

#ifndef CPP_FUNCTIONS_HPP
#define CPP_FUNCTIONS_HPP

#include "Globals.hpp"

#include <cmath>
#include <iostream>

using namespace std;

double interpol(double &x, double &x0, double &x1, double &y0, double &y1);

double npp(double _ma, double t, bool deriv);

double computeHrbmNpp(double _ma, double t, bool deriv);

double compute4boxNpp(double _ma, bool deriv);

double myLog(const double &value);

double fasC(double &ca, double &dp_co2s);

double fasT(double &rf_tot, double &t);

double rfCo2(double &_ma);

double rfEqC02Ma(double rf_co2);

double computeDpCo2s(double d_dic, double &t, bool deriv);

#endif //CPP_FUNCTIONS_HPP
