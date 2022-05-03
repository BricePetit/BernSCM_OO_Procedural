"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
import copy as cp
import sys

import numpy as np

import Globals
import functions


def computePropTempDependent(pp, s, t, q):
    """
    Wrapper that updates propagators for T-dependent IRF coefficients.

    :param pp:  IRF coefficients.
    :type pp:   PIRF.
    :param s:   IRF sensitivities.
    :type s:    SIRF.
    :param t:   Temperature perturbation.
    :type t:    Float.
    :param q:   Propagators to be updated.
    :type q:    Prop.
    """
    p = cp.deepcopy(pp)
    if t > s.t_max and not Globals.exceed:
        print("warning: temperature parametrization range of IRF ", p.name, " exceeded (",
              s.t_max, "K).", file=sys.stdout)
        Globals.exceed = True
    p.weight[:p.nscale + 1] = p.weight[:p.nscale + 1] * np.exp(s.weight[:p.nscale + 1] * t, dtype=np.longdouble)
    p.weight[:p.nscale + 1] = p.weight[:p.nscale + 1] / np.sum(p.weight[:p.nscale + 1], dtype=np.longdouble)

    p.t_scale[:p.nscale] = p.t_scale[:p.nscale] * np.exp(-s.t_scale[:p.nscale] * t, dtype=np.longdouble)

    computePropagators(p, q)


def computePropagators(p, q):
    """
    Calculate propagators.

    :param p:   IRF coefficients.
    :type p:    Pirf.
    :param q:   Propagators to be updated.
    :type q:    Prop.
    """
    if not Globals.LINEAR:
        for i in range(p.nscale):
            if p.t_scale[i] <= Globals.TEQUIL:
                q.propf[i] = p.weight[i] * p.t_scale[i]
                q.propm[i] = 0e0
            else:
                q.propf[i] = p.weight[i] * p.t_scale[i] \
                             * (1e0 - np.exp(-Globals.DELTA_T / p.t_scale[i], dtype=np.longdouble))
                q.propm[i] = np.exp(-Globals.DELTA_T / p.t_scale[i], dtype=np.longdouble)
        q.propf[p.nscale] = Globals.DELTA_T * p.weight[p.nscale]
        q.propm[p.nscale] = 1e0
    elif Globals.LINEAR:
        for i in range(p.nscale):
            q.propf[i] = p.weight[i] * p.t_scale[i] \
                         * (Globals.DELTA_T - p.t_scale[i]
                            * (1e0 - np.exp(-Globals.DELTA_T / p.t_scale[i], dtype=np.longdouble))) / Globals.DELTA_T
            q.propfo[i] = p.weight[i] * p.t_scale[i] * (
                    1e0 - np.exp(-Globals.DELTA_T / p.t_scale[i], dtype=np.longdouble)) - q.propf[i]
        q.propf[p.nscale] = p.weight[p.nscale] * Globals.DELTA_T / 2e0
        q.propfo[p.nscale] = p.weight[p.nscale] * Globals.DELTA_T - q.propf[p.nscale]

        q.propm[:p.nscale] = np.exp(- Globals.DELTA_T / p.t_scale[:p.nscale], dtype=np.longdouble)
        q.propm[p.nscale] = 1e0
    q.nscale = p.nscale


def setForcing(n):
    """
    Set forcing value and verify budget cases.

    :param n:   Time counter.
    :type n:    Int.
    """
    Globals.temp[n] = functions.interpol(
        Globals.time[n], Globals.forcing[Globals.itime[n], Globals.JTIME],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTIME],
        Globals.forcing[Globals.itime[n], Globals.JTEMP],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTEMP]
    )
    Globals.rfnc[n] = functions.interpol(
        Globals.time[n], Globals.forcing[Globals.itime[n], Globals.JTIME],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTIME],
        Globals.forcing[Globals.itime[n], Globals.JRFNC],
        Globals.forcing[Globals.itime[n] + 1, Globals.JRFNC]
    )
    Globals.rfb[n] = functions.interpol(
        Globals.time[n], Globals.forcing[Globals.itime[n], Globals.JTIME],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTIME],
        Globals.forcing[Globals.itime[n], Globals.JRFB],
        Globals.forcing[Globals.itime[n] + 1, Globals.JRFB]
    )
    Globals.ma[n] = functions.interpol(
        Globals.time[n], Globals.forcing[Globals.itime[n], Globals.JTIME],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTIME],
        Globals.forcing[Globals.itime[n], Globals.JACO2],
        Globals.forcing[Globals.itime[n] + 1, Globals.JACO2]
    )
    Globals.e_co2[n] = functions.interpol(
        Globals.time[n], Globals.forcing[Globals.itime[n], Globals.JTIME],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTIME],
        Globals.forcing[Globals.itime[n], Globals.JECO2],
        Globals.forcing[Globals.itime[n] + 1, Globals.JECO2]
    )
    Globals.fb[n] = functions.interpol(
        Globals.time[n], Globals.forcing[Globals.itime[n], Globals.JTIME],
        Globals.forcing[Globals.itime[n] + 1, Globals.JTIME],
        Globals.forcing[Globals.itime[n], Globals.JFB],
        Globals.forcing[Globals.itime[n] + 1, Globals.JFB]
    )
    if not (Globals.ma[n] == Globals.NA) or (Globals.fb[n] == Globals.NA):
        Globals.f_budget = True
    else: 
        Globals.f_budget = False
    if not (Globals.temp[n] == Globals.NA) or (Globals.rfb[n] == Globals.NA):
        Globals.rf_budget = True
    else:
        Globals.rf_budget = False

    # Check budget cases (setting RF_CO2 is not implemented)
    if (Globals.e_co2[n] == Globals.NA) or (Globals.rfnc[n] == Globals.NA):
        raise Exception(
            "eCO2 and RF_nonC must always be set, use budget_RF and budget_sink to solve for RF/emissions {}".format(
                Globals.time[n]))

    Globals.co2_budget = False
    if Globals.rf_budget and not Globals.f_budget:
        if Globals.temp[n] == Globals.NA:
            raise Exception("glob_temp_dev not set when solving for budget_RF at year {}".format(Globals.time[n]))
    elif not Globals.rf_budget and Globals.f_budget:
        if Globals.ma[n] == Globals.NA:
            raise Exception("CO2_atm not set when solving for budget_sink at year {}".format(Globals.time[n]))
    elif Globals.rf_budget and Globals.f_budget:
        if Globals.temp[n] == Globals.NA:
            raise Exception("glob_temp_dev not set when solving for budget_RF at year {}".format(Globals.time[n]))
        if Globals.ma[n] == Globals.NA:
            Globals.co2_budget = True


def stepPulse(n, f, mk, q, x):
    """
    Function to advance the integration of a tracer.

    :param n:   Time index.
    :type n:    Int.
    :param f:   Flux to mixed layer.
    :type f:    Numpy array of float.
    :param mk:  Boxes/tracer pools (input/output).
    :type mk:   Numpy array of float.
    :param q:   Propagators.
    :type q:    Prop.
    :param x:   Variable-specific multiplier/conversion factor.
    :type x:    Float.

    :return:    Return the committed temperature change (K).
    :rtype:     Float.
    """
    for j in range(q.nscale + 1):
        mk[j] = mk[j] * q.propm[j] + f[n] * q.propf[j] * x
        if Globals.LINEAR:
            mk[j] = mk[j] + f[n - 1] * q.propfo[j] * x
    return np.sum(mk[:q.nscale + 1], dtype=np.longdouble)


def timeStep(n):
    """
    Progress in the simulation step by step.

    :param n:   Correspond to the number of the current step.
    :type n:    Int.
    """
    nenner_u = 0e0
    nenner_v = 0e0
    nenner_w = 0e0
    df_npp_dma = 0e0
    ma_eq = 0e0
    dm_ao = 0e0
    e = 0e0
    if Globals.IMPLICIT_O:
        Globals.fh[n] = Globals.fh[n - 1]
        t_com = stepPulse(n, Globals.fh, Globals.tempk, Globals.ocean_prop, Globals.om_t)
    else:
        t_com = stepPulse(n - 1, Globals.fh, Globals.tempk, Globals.ocean_prop, Globals.om_t)

    # Ocean heat uptake
    if Globals.rf_budget:
        # Update heat uptake (const. flux commitment) 
        Globals.fh[n] = Globals.fh[n - 1] + (Globals.temp[n] - t_com) / (
                Globals.om_t * np.sum(Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1], dtype=np.longdouble))
        # Update Tempk!
        Globals.tempk[:Globals.ocean_prop.nscale + 1] = Globals.tempk[:Globals.ocean_prop.nscale + 1] \
                                                         + (Globals.fh[n] - Globals.fh[n - 1]) \
                                                         * Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1] \
                                                         * Globals.om_t
        # Current RF (W/m²)
        Globals.rf[n] = Globals.fh[n] / (Globals.aoc / Globals.OFRAC / Globals.PETA) \
                        + Globals.temp[n] / Globals.t2x * Globals.RF2X

        if Globals.co2_budget:
            # Calculate equivalent atmospheric CO2 (in GtC) from RF
            # It gives the atmospheric CO₂ (Gt).
            Globals.ma[n] = functions.rfEqC02Ma(Globals.rf[n] - Globals.rfnc[n])

    if Globals.LAND_T_DEP:
        if Globals.t_dep:
            if Globals.rf_budget:
                computePropTempDependent(Globals.land_pirf, Globals.land_sirf,
                                         (Globals.temp[n - 1] + Globals.temp[n]) / 2e0, Globals.land_prop)
            else:
                computePropTempDependent(Globals.land_pirf, Globals.land_sirf, (Globals.temp[n - 1] + t_com) / 2e0,
                                         Globals.land_prop)

    # land C exchange
    if Globals.f_budget:
        # solve for net C emissions
        if Globals.LINEAR:
            # (endyear value)
            if Globals.rf_budget:
                # Use actual T
                Globals.fnpp[n] = functions.npp(Globals.ma[n], Globals.temp[n], False)
            else:
                # Use committed T
                Globals.fnpp[n] = functions.npp(Globals.ma[n], t_com, False)
        else:
            # (Midyear value)
            if Globals.rf_budget:
                Globals.fnpp[n] = functions.npp((Globals.ma[n] + Globals.ma[n - 1]) / 2e0,
                                                (Globals.temp[n] + Globals.temp[n - 1]) / 2e0, False)
            else:
                Globals.fnpp[n] = functions.npp((Globals.ma[n] + Globals.ma[n - 1]) / 2e0, Globals.temp[n - 1], False)
    else:
        # Anthro emissions
        e = (Globals.e_co2[n] + Globals.e_co2[n - 1]) / 2e0 - (Globals.fb[n] + Globals.fb[n - 1]) / 2e0
        if Globals.IMPLICIT_L:
            # Auxiliary parameters
            # Commitment step with previous flux
            Globals.fnpp[n] = Globals.fnpp[n - 1]
            df_npp_dma = functions.npp(Globals.ma[n - 1], Globals.temp[n - 1], True)
            nenner_v = (df_npp_dma * np.sum(Globals.land_prop.propf[:Globals.land_prop.nscale + 1],
                                            dtype=np.longdouble) + 1e0)

    if Globals.IMPLICIT_L:
        ml_com = stepPulse(n, Globals.fnpp, Globals.mlk, Globals.land_prop, 1e0)
    else:
        if Globals.f_budget:
            ml_com = stepPulse(n, Globals.fnpp, Globals.mlk, Globals.land_prop, 1e0)
        else:
            ml_com = stepPulse(n - 1, Globals.fnpp, Globals.mlk, Globals.land_prop, 1e0)


    if Globals.IMPLICIT_O:
        # Commitment step with current flux=0
        Globals.fo[n] = 0
        dm_ao = functions.computeDpCo2s(Globals.ms[n - 1] * Globals.om_c, t_com, True) * Globals.PPMTOGT * Globals.om_c
        ma_eq = (functions.computeDpCo2s(Globals.ms[n - 1] * Globals.om_c, t_com,
                                         False) + Globals.co2_atm0) * Globals.PPMTOGT
        nenner_u = (Globals.kg_aoc * dm_ao
                    * np.sum(Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1], dtype=np.longdouble) + 1e0)
        nenner_w = Globals.DELTA_T * Globals.kg_aoc
    else:
        Globals.fo[n] = functions.fasC(Globals.ma[n - 1], Globals.dpcs[n - 1])

    ms_com = stepPulse(n, Globals.fo, Globals.msk, Globals.ocean_prop, 1e0)


    if Globals.IMPLICIT_L:
        if Globals.f_budget:
            Globals.ml[n] = ml_com
        else:
            # Implicit step for flux change (zeroE commitment for ocean)
            df_npp = df_npp_dma / (nenner_u * nenner_v + nenner_w) * (
                    Globals.ml[n - 1] - ml_com + Globals.DELTA_T * e
                    + Globals.DELTA_T * Globals.kg_aoc
                    * (+ma_eq - Globals.ma[n - 1] + dm_ao
                       * (ms_com - Globals.ms[n - 1]
                          + np.sum(Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1], dtype=np.longdouble)
                          * ((Globals.ml[n - 1] - ml_com) / Globals.DELTA_T + e))))
            Globals.ml[n] = df_npp * np.sum(Globals.land_prop.propf[:Globals.land_prop.nscale + 1], dtype=np.longdouble) + ml_com
    else:
        Globals.ml[n] = ml_com

    if Globals.IMPLICIT_O:
        # Implicit ocean step
        if Globals.f_budget:
            Globals.fo[n] = Globals.kg_aoc * (Globals.ma[n] - ma_eq - dm_ao * (ms_com - Globals.ms[n - 1])) / nenner_u
        else:
            Globals.fo[n] = Globals.kg_aoc / (nenner_u + nenner_w) \
                            * (Globals.ma[n - 1]
                               - ma_eq - dm_ao * (ms_com - Globals.ms[n - 1]) - (
                                       Globals.ml[n] - Globals.ml[n - 1]) + Globals.DELTA_T * e)
        Globals.msk[:Globals.ocean_prop.nscale + 1] = Globals.msk[:Globals.ocean_prop.nscale + 1] + Globals.fo[n] \
                                                       * Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1]
        Globals.ms[n] = np.sum(Globals.msk[:Globals.ocean_prop.nscale + 1], dtype=np.longdouble)
    else:
        Globals.ms[n] = ms_com

    # Total C budget
    if not Globals.f_budget:
        if Globals.LINEAR:
            Globals.ma[n] = Globals.ma[n - 1] + Globals.DELTA_T * (e - (Globals.fo[n - 1] + Globals.fo[n]) / 2e0) \
                            - (Globals.ml[n] - Globals.ml[n - 1])
        else:
            Globals.ma[n] = Globals.ma[n - 1] + Globals.DELTA_T * (e - Globals.fo[n]) \
                            - (Globals.ml[n] - Globals.ml[n - 1])
    # update CO₂ RF
    Globals.rfc[n] = functions.rfCo2(Globals.ma[n])

    if not Globals.rf_budget:
        # Update total RF (W/m²)
        Globals.rf[n] = Globals.rfc[n] + Globals.rfnc[n] + Globals.rfb[n]
        if Globals.t2x > 0e0:
            if Globals.IMPLICIT_O:
                # Const flux commitment
                Globals.fh[n] = (Globals.rf[n] - Globals.RF2X * t_com / Globals.t2x + Globals.fh[n - 1]
                                 * Globals.om_t
                                 * np.sum(Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1],
                                          dtype=np.longdouble)
                                 * Globals.RF2X / Globals.t2x) \
                                / (Globals.RF2X / Globals.t2x * Globals.om_t
                                   * np.sum(Globals.ocean_prop.propf[:Globals.ocean_prop.nscale + 1],
                                            dtype=np.longdouble)
                                   + Globals.OFRAC * Globals.PETA / Globals.aoc)
                Globals.tempk[:Globals.ocean_prop.nscale + 1] = Globals.tempk[:Globals.ocean_prop.nscale + 1] \
                                                                 + (Globals.fh[n] - Globals.fh[n - 1]) \
                                                                 * Globals.ocean_prop.propf[
                                                                   :Globals.ocean_prop.nscale + 1] * Globals.om_t

                Globals.temp[n] = np.sum(Globals.tempk[:Globals.ocean_prop.nscale + 1], dtype=np.longdouble)
            else:
                Globals.temp[n] = t_com
                Globals.fh[n] = functions.fasT(Globals.rf[n], Globals.temp[n])
        else:
            Globals.fh[n] = 0e0
            Globals.temp[n] = 0e0
            Globals.temp_k = np.zeros(Globals.tempk.size, dtype=np.longdouble)
    Globals.fnpp[n] = functions.npp(Globals.ma[n], Globals.temp[n], False)
    if not Globals.f_budget:
        # Update mLk with updated NPP
        Globals.mlk[:Globals.land_prop.nscale + 1] = Globals.mlk[:Globals.land_prop.nscale + 1] \
                                                      + (Globals.fnpp[n] - Globals.fnpp[n - 1]) \
                                                      * Globals.land_prop.propf[:Globals.land_prop.nscale + 1]
    Globals.dpcs[n] = functions.computeDpCo2s(Globals.ms[n] * Globals.om_c, Globals.temp[n], False)
