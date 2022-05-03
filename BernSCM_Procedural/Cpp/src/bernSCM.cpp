//
// Created by Brice Petit on 25-02-22.
//

#include "bernSCM.hpp"

/**
 * Wrapper that updates propagators for T-dependent IRF coefficients.
 *
 * @param pp    IRF coefficients.
 * @param s     IRF sensitivities.
 * @param t     Temperature perturbation.
 * @param q     Propagators to be updated.
 */
void computePropTempDependent(Pirf &pp, Sirf &s, double t, Prop &q) {
    int i;
    // Copy the struct
    Pirf p;
    p.name = pp.name;
    p.nscale = pp.nscale;
    p.weight = new double[WEIGHT_SIZE];
    memcpy(p.weight, pp.weight, sizeof(double) * WEIGHT_SIZE);
    p.t_scale = new double[T_SCALE_SIZE];
    memcpy(p.t_scale, pp.t_scale, sizeof(double) * T_SCALE_SIZE);

    double tmp_sum;
    if (t > s.t_max && !exceed) {
        cerr << "warning: temperature parametrization range of IRF " << p.name << " exceeded ("
             << s.t_max << "K)." << endl;
        exceed = true;
    }
    tmp_sum = 0;
    for (i = 0; i < p.nscale + 1; i++) {
        p.weight[i] = p.weight[i] * exp(s.weight[i] * t);
        tmp_sum += p.weight[i];
    }
    for (i = 0; i < p.nscale + 1; i++) {
        p.weight[i] = p.weight[i] / tmp_sum;
    }
    for (i = 0; i < p.nscale; i++) {
        p.t_scale[i] = p.t_scale[i] * exp(-s.t_scale[i] * t);
    }
    computePropagators(p, q);
}

/**
 * Calculate propagators.
 *
 * @param p IRF coefficients.
 * @param q Propagators to be updated.
 */
void computePropagators(const Pirf &p, Prop &q) {
    if (!LINEAR) {
        for (int i = 0; i < p.nscale; i++) {
            if (p.t_scale[i] <= TEQUIL) {
                q.propf[i] = p.weight[i] * p.t_scale[i];
                q.propm[i] = 0e0;
            } else {
                q.propf[i] = p.weight[i] * p.t_scale[i] * (1e0 - exp(-DELTA_T / p.t_scale[i]));
                q.propm[i] = exp(-DELTA_T / p.t_scale[i]);
            }
        }
        q.propf[p.nscale] = DELTA_T * p.weight[p.nscale];
        q.propm[p.nscale] = 1e0;
    } else if (LINEAR) {
        for (int i = 0; i < p.nscale; i++) {
            q.propf[i] = p.weight[i] * p.t_scale[i] *
                         (DELTA_T - p.t_scale[i] * (1e0 - exp(-DELTA_T / p.t_scale[i]))) / DELTA_T;
            q.propfo[i] = p.weight[i] * p.t_scale[i] * (1e0 - exp(-DELTA_T / p.t_scale[i])) - q.propf[i];
        }
        q.propf[p.nscale] = p.weight[p.nscale] * DELTA_T / 2e0;
        q.propfo[p.nscale] = p.weight[p.nscale] * DELTA_T - q.propf[p.nscale];

        for (int i = 0; i < p.nscale; i++) {
            q.propm[i] = exp(-DELTA_T / p.t_scale[i]);
        }
        q.propm[p.nscale] = 1e0;
    }
    q.nscale = p.nscale;
}

/**
 * Set forcing value and verify budget cases.
 *
 * @param n Time counter.
 */
void setForcing(int n) {
    temp[n] = interpol(_time[n], forcing[itime[n]][JTIME], forcing[itime[n] + 1][JTIME], forcing[itime[n]][JTEMP],
                       forcing[itime[n] + 1][JTEMP]);
    rfnc[n] = interpol(_time[n], forcing[itime[n]][JTIME], forcing[itime[n] + 1][JTIME], forcing[itime[n]][JRFNC],
                       forcing[itime[n] + 1][JRFNC]);
    rfb[n] = interpol(_time[n], forcing[itime[n]][JTIME], forcing[itime[n] + 1][JTIME], forcing[itime[n]][JRFB],
                      forcing[itime[n] + 1][JRFB]);
    ma[n] = interpol(_time[n], forcing[itime[n]][JTIME], forcing[itime[n] + 1][JTIME], forcing[itime[n]][JACO2],
                     forcing[itime[n] + 1][JACO2]);
    e_co2[n] = interpol(_time[n], forcing[itime[n]][JTIME], forcing[itime[n] + 1][JTIME], forcing[itime[n]][JECO2],
                        forcing[itime[n] + 1][JECO2]);
    fb[n] = interpol(_time[n], forcing[itime[n]][JTIME], forcing[itime[n] + 1][JTIME], forcing[itime[n]][JFB],
                     forcing[itime[n] + 1][JFB]);
    if ((ma[n] != NA) || (fb[n] == NA)) {
        f_budget = true;
    } else {
        f_budget = false;
    }
    if ((temp[n] != NA) || (rfb[n] == NA)) {
        rf_budget = true;
    } else {
        rf_budget = false;
    }
    // Check budget cases (setting RF_CO2 is not implemented)
    if ((e_co2[n] == NA) || (rfnc[n] == NA)) {
        throw "eCO2 and RF_nonC must always be set, use budget_RF and budget_sink to solve for RF/emissions " +
              to_string(_time[n]);
    }
    co2_budget = false;
    if (rf_budget && !f_budget) {
        if (temp[n] == NA) {
            throw "glob_temp_dev not set when solving for budget_RF at year " + to_string(_time[n]);
        }
    } else if (!rf_budget && f_budget) {
        if (ma[n] == NA) {
            throw "CO2_atm not set when solving for budget_sink at year " + to_string(_time[n]);
        }
    } else if (rf_budget && f_budget) {
        if (temp[n] == NA) {
            throw "glob_temp_dev not set when solving for budget_RF at year " + to_string(_time[n]);
        }
        if (ma[n] == NA) {
            co2_budget = true;
        }
    }
}

/**
 * Function to advance the integration of a tracer.
 *
 * @param n     Time index.
 * @param f     Flux to mixed layer.
 * @param mk    Boxes/tracer pools (input/output).
 * @param q     Propagators.
 * @param x     Variable-specific multiplier/conversion factor.
 * @return      Return the committed temperature change (K).
 */
double stepPulse(int n, const double *f, double *mk, Prop &q, double x) {
    double sum = 0;
    for (int j = 0; j < q.nscale + 1; j++) {
        mk[j] = mk[j] * q.propm[j] + f[n] * q.propf[j] * x;
        if (LINEAR) {
            mk[j] = mk[j] + f[n - 1] * q.propfo[j] * x;
        }
        sum += mk[j];
    }
    return sum;
}

/**
 * Progress in the simulation step by step.
 *
 * @param n Correspond to the number of the current step.
 */
void timeStep(int n) {
    double nenner_u = 0e0, nenner_v = 0e0, nenner_w = 0e0, df_npp_dma = 0e0, ma_eq = 0e0, dm_ao = 0e0, e = 0e0, t_com;
    double tmp_sum, ml_com, ms_com, df_npp;
    int a; //counter
    if (IMPLICIT_O) {
        fh[n] = fh[n - 1];
        t_com = stepPulse(n, fh, tempk, ocean_prop, om_t);
    } else {
        t_com = stepPulse(n - 1, fh, tempk, ocean_prop, om_t);
    }

    // Ocean heat uptake
    if (rf_budget) {
        // Update heat uptake (const. flux commitment)
        tmp_sum = 0;
        for (a = 0; a < ocean_prop.nscale + 1; a++) {
            tmp_sum += ocean_prop.propf[a];
        }
        fh[n] = fh[n - 1] + (temp[n] - t_com) / (om_t * tmp_sum);
        // Update Tempk!
        for (a = 0; a < ocean_prop.nscale + 1; a++) {
            tempk[a] += (fh[n] - fh[n - 1]) * ocean_prop.propf[a] * om_t;
        }
        // Current RF (W/m²)
        rf[n] = fh[n] / (aoc / OFRAC / PETA) + temp[n] / t2x * RF2X;
        if (co2_budget) {
            // Calculate equivalent atmospheric CO2 (in GtC) from RF
            // It gives the atmospheric CO₂ (Gt).
            ma[n] = rfEqC02Ma(rf[n] - rfnc[n]);
        }
    }
    if (LAND_T_DEP) {
        if (t_dep) {
            if (rf_budget) {
                computePropTempDependent(land_pirf, land_sirf, (temp[n - 1] + temp[n]) / 2e0, land_prop);
            } else {
                computePropTempDependent(land_pirf, land_sirf, (temp[n - 1] + t_com) / 2e0, land_prop);
            }
        }
    }
    // land C exchange
    if (f_budget) {
        // solve for net C emissions
        if (LINEAR) {
            // (endyear value)
            if (rf_budget) {
                // Use actual T
                fnpp[n] = npp(ma[n], temp[n], false);
            } else {
                // Use committed T
                fnpp[n] = npp(ma[n], t_com, false);
            }
        } else {
            // (Midyear value)
            if (rf_budget) {
                fnpp[n] = npp((ma[n] + ma[n - 1]) / 2e0, (temp[n] + temp[n - 1]) / 2e0, false);
            } else {
                fnpp[n] = npp((ma[n] + ma[n - 1]) / 2e0, temp[n - 1], false);
            }
        }
    } else {
        // Anthro emissions
        e = (e_co2[n] + e_co2[n - 1]) / 2e0 - (fb[n] + fb[n - 1]) / 2e0;
        if (IMPLICIT_L) {
            // Auxiliary parameters
            // Commitment step with previous flux
            fnpp[n] = fnpp[n - 1];
            df_npp_dma = npp(ma[n - 1], temp[n - 1], true);
            tmp_sum = 0;
            for (a = 0; a < land_prop.nscale + 1; a++) {
                tmp_sum += land_prop.propf[a];
            }
            nenner_v = (df_npp_dma * tmp_sum + 1e0);
        }
    }
    if (IMPLICIT_L) {
        ml_com = stepPulse(n, fnpp, mlk, land_prop, 1e0);
    } else {
        if (f_budget) {
            ml_com = stepPulse(n, fnpp, mlk, land_prop, 1e0);
        } else {
            ml_com = stepPulse(n - 1, fnpp, mlk, land_prop, 1e0);
        }
    }
    if (IMPLICIT_O) {
        // Commitment step with current flux=0
        fo[n] = 0;
        dm_ao = computeDpCo2s(ms[n - 1] * om_c, t_com, true) * PPMTOGT * om_c;
        ma_eq = (computeDpCo2s(ms[n - 1] * om_c, t_com, false) + co2_atm0) * PPMTOGT;
        tmp_sum = 0;
        for (a = 0; a < ocean_prop.nscale + 1; a++) {
            tmp_sum += ocean_prop.propf[a];
        }
        nenner_u = (kg_aoc * dm_ao * tmp_sum + 1e0);
        nenner_w = DELTA_T * kg_aoc;
    } else {
        fo[n] = fasC(ma[n - 1], dpcs[n - 1]);
    }

    ms_com = stepPulse(n, fo, msk, ocean_prop, 1e0);

    if (IMPLICIT_L) {
        if (f_budget) {
            ml[n] = ml_com;
        } else {
            // Implicit step for flux change (zeroE commitment for ocean)
            tmp_sum = 0;
            for (a = 0; a < ocean_prop.nscale + 1; a++) {
                tmp_sum += ocean_prop.propf[a];
            }
            df_npp = df_npp_dma / (nenner_u * nenner_v + nenner_w) *
                    (ml[n - 1] - ml_com + DELTA_T * e + DELTA_T * kg_aoc *
                    (+ma_eq - ma[n - 1] + dm_ao * (ms_com - ms[n - 1] + tmp_sum *
                    ((ml[n - 1] - ml_com) / DELTA_T + e))));
            tmp_sum = 0;
            for (a = 0; a < land_prop.nscale + 1; a++) {
                tmp_sum += land_prop.propf[a];
            }
            ml[n] = df_npp * tmp_sum + ml_com;
        }
    } else {
        ml[n] = ml_com;
    }

    if (IMPLICIT_O) {
        // Implicit ocean step
        if (f_budget) {
            fo[n] = kg_aoc * (ma[n] - ma_eq - dm_ao * (ms_com - ms[n - 1])) / nenner_u;
        } else {
            fo[n] = kg_aoc / (nenner_u + nenner_w) * (ma[n - 1] - ma_eq - dm_ao *
                    (ms_com - ms[n - 1]) - (ml[n] - ml[n - 1]) + DELTA_T * e);
        }
        tmp_sum = 0;
        for (a = 0; a < ocean_prop.nscale + 1; a++) {
            msk[a] = msk[a] + fo[n] * ocean_prop.propf[a];
            tmp_sum += msk[a];
        }
        ms[n] = tmp_sum;
    } else {
        ms[n] = ms_com;
    }

    // Total C budget
    if (!f_budget) {
        if (LINEAR) {
            ma[n] = ma[n - 1] + DELTA_T * (e - (fo[n - 1] + fo[n]) / 2e0) - (ml[n] - ml[n - 1]);
        } else {
            ma[n] = ma[n - 1] + DELTA_T * (e - fo[n]) - (ml[n] - ml[n - 1]);
        }
    }
    // update CO₂ RF
    rfc[n] = rfCo2(ma[n]);

    if (!rf_budget) {
        // Update total RF (W/m²)
        rf[n] = rfc[n] + rfnc[n] + rfb[n];
        if (t2x > 0e0) {
            if (IMPLICIT_O) {
                // Const flux commitment
                tmp_sum = 0;
                for (a = 0; a < ocean_prop.nscale + 1; a++) {
                    tmp_sum += ocean_prop.propf[a];
                }
                fh[n] = (rf[n] - RF2X * t_com / t2x + fh[n - 1] * om_t * tmp_sum * RF2X / t2x) /
                        (RF2X / t2x * om_t * tmp_sum + OFRAC * PETA / aoc);
                tmp_sum = 0;
                for (a = 0; a < ocean_prop.nscale + 1; a++) {
                    tempk[a] += (fh[n] - fh[n - 1]) * ocean_prop.propf[a] * om_t;
                    tmp_sum += tempk[a];
                }
                temp[n] = tmp_sum;
            } else {
                temp[n] = t_com;
                fh[n] = fasT(rf[n], temp[n]);
            }
        } else {
            fh[n] = 0e0;
            temp[n] = 0e0;
            delete[] tempk;
            tempk = new double[ocean_pirf.nscale + 1]{0};
        }
    }
    fnpp[n] = npp(ma[n], temp[n], false);
    if (!f_budget) {
        //Update mLk with updated NPP
        for (a = 0; a < land_prop.nscale + 1; a++) {
            mlk[a] = mlk[a] + (fnpp[n] - fnpp[n - 1]) * land_prop.propf[a];
        }
    }
    dpcs[n] = computeDpCo2s(ms[n] * om_c, temp[n], false);
}
