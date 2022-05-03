//
// Created by Brice Petit on 25-02-22.
//

#ifndef CPP_BERNSCM_HPP
#define CPP_BERNSCM_HPP

#include "functions.hpp"
#include "Globals.hpp"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace std;

using namespace std;

void computePropTempDependent(Pirf &pp, Sirf &s, double t, Prop &q);

void computePropagators(const Pirf& p, Prop &q);

void setForcing(int n);

double stepPulse(int n, const double *f, double *mk, Prop &q, double x);

void timeStep(int n);

#endif //CPP_BERNSCM_HPP
