/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "bernSCM_init.hpp"

/**
 * Initialization for the simulation.
 */
void initialize() {
    int nin = FORCING_ROWS - 1, j = 0;
    double tin[FORCING_ROWS];
    if (LINEAR && TEQUIL > 0) {
        throw "Equilibrated time steps not implemented for linear.";
    } else if (LINEAR && !(IMPLICIT_O || IMPLICIT_L) == 1) {
        throw "Linear scheme only implicit.";
    }
    for (int i = 0; i < FORCING_ROWS; i++) {
        tin[i] = forcing[i][JTIME];
    }
    ntime = int((tin[nin] - tin[0]) / DELTA_T + 1e-6) + 1;
    _time = new double[ntime]{0};
    itime = new long int[ntime]{0};

    for (int i = 0; i < ntime; i++) {
        _time[i] = tin[0] + i * DELTA_T;
        while (_time[i] > tin[j + 1]) {
            j += 1;
        }
        itime[i] = j;
    }

    rf = new double[ntime]{0};
    rfnc = new double[ntime]{0};
    rfc = new double[ntime]{0};
    rfb = new double[ntime]{0};
    fh = new double[ntime]{0};
    fb = new double[ntime]{0};
    fo = new double[ntime]{0};
    ma = new double[ntime]{0};
    e_co2 = new double[ntime]{0};
    fnpp = new double[ntime]{NPP0};
    dpcs = new double[ntime]{0};

    computePropagators(land_pirf, land_prop);
    computePropagators(ocean_pirf, ocean_prop);

    // Initialize the stock of carbon for the land.
    mlk = new double[land_pirf.nscale + 1]{0};
    ml = new double[ntime]{0};
    for (int i = 0; i < land_pirf.nscale; i++) {
        mlk[i] = (NPP0 * land_pirf.weight[i] * land_pirf.t_scale[i]);
        ml[0] += mlk[i];
    }

    // Initial ocean mixed layer C stock perturbation
    msk = new double[ocean_pirf.nscale + 1]{0};
    ms = new double[ntime]{0};

    // Initial ocean temperature perturbation
    tempk = new double[ocean_pirf.nscale + 1]{0};
    temp = new double[ntime]{0};
}